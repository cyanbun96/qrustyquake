// Copyright (C) 1996-2001 Id Software, Inc.
// Copyright (C) 2010-2014 QuakeSpasm developers
// GPLv3 See LICENSE for details.

#include "quakedef.h"

#define sfunc net_drivers[sock->driver] // readability macros
#define dfunc net_drivers[net_driverlevel]

s32 net_driverlevel;
f64 net_time;
s32 hostCacheCount = 0;
hostcache_t hostcache[HOSTCACHESIZE];
qsocket_t *net_activeSockets = NULL;
qsocket_t *net_freeSockets = NULL;
s32 net_numsockets = 0;
bool tcpipAvailable = 0;
s32 net_hostport;
s32 DEFAULTnet_hostport = 26000;
s8 my_ipx_address[NET_NAMELEN];
s8 my_tcpip_address[NET_NAMELEN];
static bool listening = 0;
bool slistInProgress = 0;
bool slistSilent = 0;
bool slistLocal = 1;
static f64 slistStartTime;
static s32 slistLastShown;
sizebuf_t net_message;
s32 net_activeconnections = 0;
s32 messagesSent = 0;
s32 messagesReceived = 0;
s32 unreliableMessagesSent = 0;
s32 unreliableMessagesReceived = 0;
static PollProcedure *pollProcedureList = NULL;
static void Slist_Send(void *);
static void Slist_Poll(void *);
static PollProcedure slistSendProcedure = { NULL, 0.0, Slist_Send, NULL };
static PollProcedure slistPollProcedure = { NULL, 0.0, Slist_Poll, NULL };

f64 SetNetTime()
{
	net_time = Sys_DoubleTime();
	return net_time;
}

// Called by drivers when a new communications endpoint is required
// The sequence and buffer fields will be filled in properly
qsocket_t *NET_NewQSocket()
{
	qsocket_t *sock;
	if (net_freeSockets == NULL)
		return NULL;
	if (net_activeconnections >= svs.maxclients)
		return NULL;
	// get one from free list
	sock = net_freeSockets;
	net_freeSockets = sock->next;
	// add it to active list
	sock->next = net_activeSockets;
	net_activeSockets = sock;
	sock->disconnected = 0;
	sock->connecttime = net_time;
	Q_strcpy(sock->address, "UNSET ADDRESS");
	sock->driver = net_driverlevel;
	sock->socket = 0;
	sock->driverdata = NULL;
	sock->canSend = 1;
	sock->sendNext = 0;
	sock->lastMessageTime = net_time;
	sock->ackSequence = 0;
	sock->sendSequence = 0;
	sock->unreliableSendSequence = 0;
	sock->sendMessageLength = 0;
	sock->receiveSequence = 0;
	sock->unreliableReceiveSequence = 0;
	sock->receiveMessageLength = 0;
	return sock;
}

void NET_FreeQSocket(qsocket_t *sock)
{
	qsocket_t *s;
	// remove it from active list
	if (sock == net_activeSockets)
		net_activeSockets = net_activeSockets->next;
	else {
		for (s = net_activeSockets; s; s = s->next) {
			if (s->next == sock) {
				s->next = sock->next;
				break;
			}
		}
		if (!s)
			Sys_Error("NET_FreeQSocket: not active");
	}
	// add it to free list
	sock->next = net_freeSockets;
	net_freeSockets = sock;
	sock->disconnected = 1;
}

static void NET_Listen_f()
{
	if (Cmd_Argc() != 2) {
		Con_Printf("\"listen\" is \"%d\"\n", listening ? 1 : 0);
		return;
	}
	listening = Q_atoi(Cmd_Argv(1)) ? 1 : 0;
	for (net_driverlevel = 0; net_driverlevel < net_numdrivers;
	     net_driverlevel++) {
		if (net_drivers[net_driverlevel].initialized == 0)
			continue;
		dfunc.Listen(listening);
	}
}

static void MaxPlayers_f()
{
	s32 n;
	if (Cmd_Argc() != 2) {
		Con_Printf("\"maxplayers\" is \"%d\"\n", svs.maxclients);
		return;
	}
	if (sv.active) {
		Con_Printf
		    ("maxplayers can not be changed while a server is running.\n");
		return;
	}
	n = Q_atoi(Cmd_Argv(1));
	if (n < 1)
		n = 1;
	if (n > svs.maxclientslimit) {
		n = svs.maxclientslimit;
		Con_Printf("\"maxplayers\" set to \"%d\"\n", n);
	}
	if ((n == 1) && listening)
		Cbuf_AddText("listen 0\n");
	if ((n > 1) && (!listening))
		Cbuf_AddText("listen 1\n");
	svs.maxclients = n;
	if (n == 1)
		Cvar_Set("deathmatch", "0");
	else
		Cvar_Set("deathmatch", "1");
}

static void NET_Port_f()
{
	s32 n;
	if (Cmd_Argc() != 2) {
		Con_Printf("\"port\" is \"%d\"\n", net_hostport);
		return;
	}
	n = Q_atoi(Cmd_Argv(1));
	if (n < 1 || n > 65534) {
		Con_Printf("Bad value, must be between 1 and 65534\n");
		return;
	}
	DEFAULTnet_hostport = n;
	net_hostport = n;
	if (listening) {
		// force a change to the new port
		Cbuf_AddText("listen 0\n");
		Cbuf_AddText("listen 1\n");
	}
}

static void PrintSlistHeader()
{
	Con_Printf("Server          Map             Users\n");
	Con_Printf("--------------- --------------- -----\n");
	slistLastShown = 0;
}

static void PrintSlist()
{
	s32 n;
	for (n = slistLastShown; n < hostCacheCount; n++) {
		if (hostcache[n].maxusers)
			Con_Printf("%-15.15s %-15.15s %2u/%2u\n",
				   hostcache[n].name, hostcache[n].map,
				   hostcache[n].users, hostcache[n].maxusers);
		else
			Con_Printf("%-15.15s %-15.15s\n", hostcache[n].name,
				   hostcache[n].map);
	}
	slistLastShown = n;
}

static void PrintSlistTrailer()
{
	if (hostCacheCount)
		Con_Printf("== end list ==\n\n");
	else
		Con_Printf("No Quake servers found.\n\n");
}

void NET_Slist_f()
{
	if (slistInProgress)
		return;
	if (!slistSilent) {
		Con_Printf("Looking for Quake servers...\n");
		PrintSlistHeader();
	}
	slistInProgress = 1;
	slistStartTime = Sys_DoubleTime();
	SchedulePollProcedure(&slistSendProcedure, 0.0);
	SchedulePollProcedure(&slistPollProcedure, 0.1);
	hostCacheCount = 0;
}

static void Slist_Send(void *unused)
{
	(void)unused;
	for (net_driverlevel = 0; net_driverlevel < net_numdrivers;
	     net_driverlevel++) {
		if (!slistLocal && IS_LOOP_DRIVER(net_driverlevel))
			continue;
		if (net_drivers[net_driverlevel].initialized == 0)
			continue;
		dfunc.SearchForHosts(1);
	}
	if ((Sys_DoubleTime() - slistStartTime) < 0.5)
		SchedulePollProcedure(&slistSendProcedure, 0.75);
}

static void Slist_Poll(void *unused)
{
	(void)unused;
	for (net_driverlevel = 0; net_driverlevel < net_numdrivers;
	     net_driverlevel++) {
		if (!slistLocal && IS_LOOP_DRIVER(net_driverlevel))
			continue;
		if (net_drivers[net_driverlevel].initialized == 0)
			continue;
		dfunc.SearchForHosts(0);
	}
	if (!slistSilent)
		PrintSlist();
	if ((Sys_DoubleTime() - slistStartTime) < 1.5) {
		SchedulePollProcedure(&slistPollProcedure, 0.1);
		return;
	}
	if (!slistSilent)
		PrintSlistTrailer();
	slistInProgress = 0;
	slistSilent = 0;
	slistLocal = 1;
}

qsocket_t *NET_Connect(const s8 *host)
{
	qsocket_t *ret;
	s32 n;
	s32 numdrivers = net_numdrivers;
	SetNetTime();
	if (host && *host == 0)
		host = NULL;
	if (host) {
		if (q_strcasecmp(host, "local") == 0) {
			numdrivers = 1;
			goto JustDoIt;
		}
		if (hostCacheCount) {
			for (n = 0; n < hostCacheCount; n++)
				if (q_strcasecmp(host, hostcache[n].name) == 0) {
					host = hostcache[n].cname;
					break;
				}
			if (n < hostCacheCount)
				goto JustDoIt;
		}
	}
	slistSilent = host ? 1 : 0;
	NET_Slist_f();
	while (slistInProgress)
		NET_Poll();
	if (host == NULL) {
		if (hostCacheCount != 1)
			return NULL;
		host = hostcache[0].cname;
		Con_Printf("Connecting to...\n%s @ %s\n\n", hostcache[0].name,
			   host);
	}
	if (hostCacheCount) {
		for (n = 0; n < hostCacheCount; n++) {
			if (q_strcasecmp(host, hostcache[n].name) == 0) {
				host = hostcache[n].cname;
				break;
			}
		}
	}
JustDoIt:
	for (net_driverlevel = 0; net_driverlevel < numdrivers;
	     net_driverlevel++) {
		if (net_drivers[net_driverlevel].initialized == 0)
			continue;
		ret = dfunc.Connect(host);
		if (ret)
			return ret;
	}
	if (host) {
		Con_Printf("\n");
		PrintSlistHeader();
		PrintSlist();
		PrintSlistTrailer();
	}
	return NULL;
}

qsocket_t *NET_CheckNewConnections()
{
	qsocket_t *ret;
	SetNetTime();
	for (net_driverlevel = 0; net_driverlevel < net_numdrivers;
	     net_driverlevel++) {
		if (net_drivers[net_driverlevel].initialized == 0)
			continue;
		if (!IS_LOOP_DRIVER(net_driverlevel) && listening == 0)
			continue;
		ret = dfunc.CheckNewConnections();
		if (ret) {
			return ret;
		}
	}
	return NULL;
}

void NET_Close(qsocket_t *sock)
{
	if (!sock)
		return;
	if (sock->disconnected)
		return;
	SetNetTime();
	// call the driver_Close function
	sfunc.Close(sock);
	NET_FreeQSocket(sock);
}

// If there is a complete message, return it in net_message
// returns 0 if no data is waiting
// returns 1 if a message was received
// returns -1 if connection is invalid
s32 NET_GetMessage(qsocket_t *sock)
{
	s32 ret;
	if (!sock)
		return -1;
	if (sock->disconnected) {
		Con_Printf("NET_GetMessage: disconnected socket\n");
		return -1;
	}
	SetNetTime();
	ret = sfunc.QGetMessage(sock);
	// see if this connection has timed out
	if (ret == 0 && !IS_LOOP_DRIVER(sock->driver)) {
		if (net_time - sock->lastMessageTime > net_messagetimeout.value) {
			NET_Close(sock);
			return -1;
		}
	}
	if (ret > 0) {
		if (!IS_LOOP_DRIVER(sock->driver)) {
			sock->lastMessageTime = net_time;
			if (ret == 1)
				messagesReceived++;
			else if (ret == 2)
				unreliableMessagesReceived++;
		}
	}
	return ret;
}

// Try to send a complete length+message unit over the reliable stream.
// returns 0 if the message cannot be delivered reliably, but the connection
// 		is still considered valid
// returns 1 if the message was sent properly
// returns -1 if the connection died
s32 NET_SendMessage(qsocket_t *sock, sizebuf_t *data)
{
	s32 r;
	if (!sock)
		return -1;
	if (sock->disconnected) {
		Con_Printf("NET_SendMessage: disconnected socket\n");
		return -1;
	}
	SetNetTime();
	r = sfunc.QSendMessage(sock, data);
	if (r == 1 && !IS_LOOP_DRIVER(sock->driver))
		messagesSent++;
	return r;
}

s32 NET_SendUnreliableMessage(qsocket_t *sock, sizebuf_t *data)
{
	s32 r;
	if (!sock)
		return -1;
	if (sock->disconnected) {
		Con_Printf("NET_SendMessage: disconnected socket\n");
		return -1;
	}
	SetNetTime();
	r = sfunc.SendUnreliableMessage(sock, data);
	if (r == 1 && !IS_LOOP_DRIVER(sock->driver))
		unreliableMessagesSent++;
	return r;
}

// Returns 1 or 0 if the given qsocket can currently accept a
// message to be transmitted.
bool NET_CanSendMessage(qsocket_t *sock)
{
	if (!sock)
		return 0;
	if (sock->disconnected)
		return 0;
	SetNetTime();
	return sfunc.CanSendMessage(sock);
}

s32 NET_SendToAll(sizebuf_t *data, f64 blocktime)
{
	f64 start;
	s32 i = 0;
	s32 count = 0;
	bool msg_init[MAX_SCOREBOARD];	/* did we write the message to the client's connection  */
	bool msg_sent[MAX_SCOREBOARD];	/* did the msg arrive its destination (canSend state).  */
	for (host_client = svs.clients; i < svs.maxclients; i++, host_client++){
		if (host_client->netconnection && host_client->active) {
			if (IS_LOOP_DRIVER(host_client->netconnection->driver)){
				NET_SendMessage(host_client->netconnection,
						data);
				msg_init[i] = 1;
				msg_sent[i] = 1;
				continue;
			}
			count++;
			msg_init[i] = 0;
			msg_sent[i] = 0;
		} else {
			msg_init[i] = 1;
			msg_sent[i] = 1;
		}
	}
	start = Sys_DoubleTime();
	while (count) {
		count = 0;
		for (i = 0, host_client = svs.clients; i < svs.maxclients;
		     i++, host_client++) {
			if (!msg_init[i]) {
				if (NET_CanSendMessage
				    (host_client->netconnection)) {
					msg_init[i] = 1;
					NET_SendMessage(host_client->
							netconnection, data);
				} else {
					NET_GetMessage(host_client->
						       netconnection);
				}
				count++;
				continue;
			}
			if (!msg_sent[i]) {
				if (NET_CanSendMessage
				    (host_client->netconnection)) {
					msg_sent[i] = 1;
				} else {
					NET_GetMessage(host_client->
						       netconnection);
				}
				count++;
				continue;
			}
		}
		if ((Sys_DoubleTime() - start) > blocktime)
			break;
	}
	return count;
}

void NET_Init()
{
	s32 i;
	qsocket_t *s;
	i = COM_CheckParm("-port");
	if (!i)
		i = COM_CheckParm("-udpport");
	if (!i)
		i = COM_CheckParm("-ipxport");
	if (i) {
		if (i < com_argc - 1)
			DEFAULTnet_hostport = Q_atoi(com_argv[i + 1]);
		else
			Sys_Error
			    ("NET_Init: you must specify a number after -port");
	}
	net_hostport = DEFAULTnet_hostport;
	net_numsockets = svs.maxclientslimit;
	if (cls.state != ca_dedicated)
		net_numsockets++;
	if (COM_CheckParm("-listen") || cls.state == ca_dedicated)
		listening = 1;
	SetNetTime();
	for (i = 0; i < net_numsockets; i++) {
		s = (qsocket_t *) Hunk_AllocName(sizeof(qsocket_t), "qsocket");
		s->next = net_freeSockets;
		net_freeSockets = s;
		s->disconnected = 1;
	}
	// allocate space for network message buffer
	SZ_Alloc(&net_message, NET_MAXMESSAGE);
	Cvar_RegisterVariable(&net_messagetimeout);
	Cvar_RegisterVariable(&hostname);
	Cmd_AddCommand("slist", NET_Slist_f);
	Cmd_AddCommand("listen", NET_Listen_f);
	Cmd_AddCommand("maxplayers", MaxPlayers_f);
	Cmd_AddCommand("port", NET_Port_f);
	// initialize all the drivers
	for (i = net_driverlevel = 0; net_driverlevel < net_numdrivers;
	     net_driverlevel++) {
		if (net_drivers[net_driverlevel].Init() == -1)
			continue;
		i++;
		net_drivers[net_driverlevel].initialized = 1;
		if (listening)
			net_drivers[net_driverlevel].Listen(1);
	}
	/* Loop_Init() returns -1 for dedicated server case,
	 * therefore the i == 0 check is correct */
	if (i == 0 && cls.state == ca_dedicated) {
		Sys_Error("Network not available!");
	}
	if (*my_ipx_address) {
		Con_DPrintf("IPX address %s\n", my_ipx_address);
	}
	if (*my_tcpip_address) {
		Con_DPrintf("TCP/IP address %s\n", my_tcpip_address);
	}
}

void NET_Shutdown()
{
	qsocket_t *sock;
	SetNetTime();
	for (sock = net_activeSockets; sock; sock = sock->next)
		NET_Close(sock);
	// shutdown the drivers
	for (net_driverlevel = 0; net_driverlevel < net_numdrivers;
	     net_driverlevel++) {
		if (net_drivers[net_driverlevel].initialized == 1) {
			net_drivers[net_driverlevel].Shutdown();
			net_drivers[net_driverlevel].initialized = 0;
		}
	}
}

void NET_Poll()
{
	PollProcedure *pp;
	SetNetTime();
	for (pp = pollProcedureList; pp; pp = pp->next) {
		if (pp->nextTime > net_time)
			break;
		pollProcedureList = pp->next;
		pp->procedure(pp->arg);
	}
}

void SchedulePollProcedure(PollProcedure *proc, f64 timeOffset)
{
	PollProcedure *pp, *prev;
	proc->nextTime = Sys_DoubleTime() + timeOffset;
	for (pp = pollProcedureList, prev = NULL; pp; pp = pp->next) {
		if (pp->nextTime >= proc->nextTime)
			break;
		prev = pp;
	}
	if (prev == NULL) {
		proc->next = pollProcedureList;
		pollProcedureList = proc;
		return;
	}
	proc->next = pp;
	prev->next = proc;
}
#undef sfunc
#undef dfunc
