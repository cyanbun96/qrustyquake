// Copyright (C) 1996-2001 Id Software, Inc.
// Copyright (C) 2010-2014 QuakeSpasm developers
// GPLv3 See LICENSE for details.
#include "quakedef.h"

#define sfunc net_landrivers[sock->landriver] // readability macros
#define dfunc net_landrivers[net_landriverlevel]

static struct {
	u32 length;
	u32 sequence;
	u8 data[MAX_DATAGRAM];
} packetBuffer;

static s32 net_landriverlevel;
/* statistic counters */
static s32 packetsSent = 0;
static s32 packetsReSent = 0;
static s32 packetsReceived = 0;
static s32 receivedDuplicateCount = 0;
static s32 shortPacketCount = 0;
static s32 droppedDatagrams;
static s32 myDriverLevel;
extern bool m_return_onerror;
extern s8 m_return_reason[32];

static s8 *StrAddr(struct qsockaddr *addr)
{
	static s8 buf[34];
	u8 *p = (u8 *) addr;
	s32 n;
	for (n = 0; n < 16; n++)
		sprintf(buf + n * 2, "%02x", *p++);
	return buf;
}

s32 Datagram_SendMessage(qsocket_t *sock, sizebuf_t *data)
{
	u32 packetLen;
	u32 dataLen;
	u32 eom;
	Q_memcpy(sock->sendMessage, data->data, data->cursize);
	sock->sendMessageLength = data->cursize;
	if (data->cursize <= MAX_DATAGRAM) {
		dataLen = data->cursize;
		eom = NETFLAG_EOM;
	} else {
		dataLen = MAX_DATAGRAM;
		eom = 0;
	}
	packetLen = NET_HEADERSIZE + dataLen;
	packetBuffer.length = BigLong(packetLen | (NETFLAG_DATA | eom));
	packetBuffer.sequence = BigLong(sock->sendSequence++);
	Q_memcpy(packetBuffer.data, sock->sendMessage, dataLen);
	sock->canSend = 0;
	if (sfunc.
	    Write(sock->socket, (u8 *) & packetBuffer, packetLen,
		  &sock->addr) == -1)
		return -1;
	sock->lastSendTime = net_time;
	packetsSent++;
	return 1;
}

static s32 SendMessageNext(qsocket_t *sock)
{
	u32 packetLen;
	u32 dataLen;
	u32 eom;
	if (sock->sendMessageLength <= MAX_DATAGRAM) {
		dataLen = sock->sendMessageLength;
		eom = NETFLAG_EOM;
	} else {
		dataLen = MAX_DATAGRAM;
		eom = 0;
	}
	packetLen = NET_HEADERSIZE + dataLen;
	packetBuffer.length = BigLong(packetLen | (NETFLAG_DATA | eom));
	packetBuffer.sequence = BigLong(sock->sendSequence++);
	Q_memcpy(packetBuffer.data, sock->sendMessage, dataLen);
	sock->sendNext = 0;
	if (sfunc.
	    Write(sock->socket, (u8 *) & packetBuffer, packetLen,
		  &sock->addr) == -1)
		return -1;
	sock->lastSendTime = net_time;
	packetsSent++;
	return 1;
}

static s32 ReSendMessage(qsocket_t *sock)
{
	u32 packetLen;
	u32 dataLen;
	u32 eom;
	if (sock->sendMessageLength <= MAX_DATAGRAM) {
		dataLen = sock->sendMessageLength;
		eom = NETFLAG_EOM;
	} else {
		dataLen = MAX_DATAGRAM;
		eom = 0;
	}
	packetLen = NET_HEADERSIZE + dataLen;
	packetBuffer.length = BigLong(packetLen | (NETFLAG_DATA | eom));
	packetBuffer.sequence = BigLong(sock->sendSequence - 1);
	Q_memcpy(packetBuffer.data, sock->sendMessage, dataLen);
	sock->sendNext = 0;
	if (sfunc.
	    Write(sock->socket, (u8 *) & packetBuffer, packetLen,
		  &sock->addr) == -1)
		return -1;
	sock->lastSendTime = net_time;
	packetsReSent++;
	return 1;
}

bool Datagram_CanSendMessage(qsocket_t *sock)
{
	if (sock->sendNext)
		SendMessageNext(sock);
	return sock->canSend;
}

bool Datagram_CanSendUnreliableMessage(qsocket_t *sock)
{
	(void)sock; // the other sock, i found it.
	return 1;
}

s32 Datagram_SendUnreliableMessage(qsocket_t *sock, sizebuf_t *data)
{
	s32 packetLen;
	packetLen = NET_HEADERSIZE + data->cursize;
	packetBuffer.length = BigLong(packetLen | NETFLAG_UNRELIABLE);
	packetBuffer.sequence = BigLong(sock->unreliableSendSequence++);
	Q_memcpy(packetBuffer.data, data->data, data->cursize);
	if (sfunc.
	    Write(sock->socket, (u8 *) & packetBuffer, packetLen,
		  &sock->addr) == -1)
		return -1;
	packetsSent++;
	return 1;
}

s32 Datagram_GetMessage(qsocket_t *sock)
{
	u32 length;
	u32 flags;
	s32 ret = 0;
	struct qsockaddr readaddr;
	u32 sequence;
	u32 count;
	if (!sock->canSend)
		if ((net_time - sock->lastSendTime) > 1.0)
			ReSendMessage(sock);
	while (1) {
		length =
		    (u32)sfunc.Read(sock->socket,
					     (u8 *) & packetBuffer,
					     NET_DATAGRAMSIZE, &readaddr);
		//      if ((rand() & 255) > 220)
		//              continue;
		if (length == 0)
			break;
		if (length == (u32)-1) {
			Con_Printf("Read error\n");
			return -1;
		}
		if (sfunc.AddrCompare(&readaddr, &sock->addr) != 0) {
			Con_Printf("Forged packet received\n");
			Con_Printf("Expected: %s\n", StrAddr(&sock->addr));
			Con_Printf("Received: %s\n", StrAddr(&readaddr));
			continue;
		}
		if (length < NET_HEADERSIZE) {
			shortPacketCount++;
			continue;
		}
		length = BigLong(packetBuffer.length);
		flags = length & (~NETFLAG_LENGTH_MASK);
		length &= NETFLAG_LENGTH_MASK;
		if (flags & NETFLAG_CTL)
			continue;
		sequence = BigLong(packetBuffer.sequence);
		packetsReceived++;
		if (flags & NETFLAG_UNRELIABLE) {
			if (sequence < sock->unreliableReceiveSequence) {
				Con_DPrintf("Got a stale datagram\n");
				ret = 0;
				break;
			}
			if (sequence != sock->unreliableReceiveSequence) {
				count =
				    sequence - sock->unreliableReceiveSequence;
				droppedDatagrams += count;
				Con_DPrintf("Dropped %u datagram(s)\n", count);
			}
			sock->unreliableReceiveSequence = sequence + 1;
			length -= NET_HEADERSIZE;
			SZ_Clear(&net_message);
			SZ_Write(&net_message, packetBuffer.data, length);
			ret = 2;
			break;
		}
		if (flags & NETFLAG_ACK) {
			if (sequence != (sock->sendSequence - 1)) {
				Con_DPrintf("Stale ACK received\n");
				continue;
			}
			if (sequence == sock->ackSequence) {
				sock->ackSequence++;
				if (sock->ackSequence != sock->sendSequence)
					Con_DPrintf("ack sequencing error\n");
			} else {
				Con_DPrintf("Duplicate ACK received\n");
				continue;
			}
			sock->sendMessageLength -= MAX_DATAGRAM;
			if (sock->sendMessageLength > 0) {
				memmove(sock->sendMessage,
					sock->sendMessage + MAX_DATAGRAM,
					sock->sendMessageLength);
				sock->sendNext = 1;
			} else {
				sock->sendMessageLength = 0;
				sock->canSend = 1;
			}
			continue;
		}
		if (flags & NETFLAG_DATA) {
			packetBuffer.length =
			    BigLong(NET_HEADERSIZE | NETFLAG_ACK);
			packetBuffer.sequence = BigLong(sequence);
			sfunc.Write(sock->socket, (u8 *) & packetBuffer,
				    NET_HEADERSIZE, &readaddr);
			if (sequence != sock->receiveSequence) {
				receivedDuplicateCount++;
				continue;
			}
			sock->receiveSequence++;
			length -= NET_HEADERSIZE;
			if (flags & NETFLAG_EOM) {
				SZ_Clear(&net_message);
				SZ_Write(&net_message, sock->receiveMessage,
					 sock->receiveMessageLength);
				SZ_Write(&net_message, packetBuffer.data,
					 length);
				sock->receiveMessageLength = 0;
				ret = 1;
				break;
			}
			Q_memcpy(sock->receiveMessage +
				 sock->receiveMessageLength, packetBuffer.data,
				 length);
			sock->receiveMessageLength += length;
			continue;
		}
	}
	if (sock->sendNext)
		SendMessageNext(sock);
	return ret;
}

static void PrintStats(qsocket_t *s)
{
	Con_Printf("canSend = %4u   \n", s->canSend);
	Con_Printf("sendSeq = %4u   ", s->sendSequence);
	Con_Printf("recvSeq = %4u   \n", s->receiveSequence);
	Con_Printf("\n");
}

static void NET_Stats_f()
{
	qsocket_t *s;
	if (Cmd_Argc() == 1) {
		Con_Printf("unreliable messages sent   = %i\n",
			   unreliableMessagesSent);
		Con_Printf("unreliable messages recv   = %i\n",
			   unreliableMessagesReceived);
		Con_Printf("reliable messages sent     = %i\n", messagesSent);
		Con_Printf("reliable messages received = %i\n",
			   messagesReceived);
		Con_Printf("packetsSent                = %i\n", packetsSent);
		Con_Printf("packetsReSent              = %i\n", packetsReSent);
		Con_Printf("packetsReceived            = %i\n",
			   packetsReceived);
		Con_Printf("receivedDuplicateCount     = %i\n",
			   receivedDuplicateCount);
		Con_Printf("shortPacketCount           = %i\n",
			   shortPacketCount);
		Con_Printf("droppedDatagrams           = %i\n",
			   droppedDatagrams);
	} else if (Q_strcmp(Cmd_Argv(1), "*") == 0) {
		for (s = net_activeSockets; s; s = s->next)
			PrintStats(s);
		for (s = net_freeSockets; s; s = s->next)
			PrintStats(s);
	} else {
		for (s = net_activeSockets; s; s = s->next) {
			if (q_strcasecmp(Cmd_Argv(1), s->address) == 0)
				break;
		}
		if (s == NULL) {
			for (s = net_freeSockets; s; s = s->next) {
				if (q_strcasecmp(Cmd_Argv(1), s->address) == 0)
					break;
			}
		}
		if (s == NULL)
			return;
		PrintStats(s);
	}
}

// recognize ip:port (based on ProQuake)
static const s8 *Strip_Port(const s8 *host)
{
	static s8 noport[MAX_QPATH];
	/* array size as in Host_Connect_f() */
	s8 *p;
	s32 port;
	if (!host || !*host)
		return host;
	Q_strncpy(noport, host, sizeof(noport));
	if ((p = Q_strrchr(noport, ':')) == NULL)
		return host;
	*p++ = '\0';
	port = Q_atoi(p);
	if (port > 0 && port < 65536 && port != net_hostport) {
		net_hostport = port;
		Con_Printf("Port set to %d\n", net_hostport);
	}
	return noport;
}

static bool testInProgress = 0;
static s32 testPollCount;
static s32 testDriver;
static sys_socket_t testSocket;
static void Test_Poll(void *);
static PollProcedure testPollProcedure = { NULL, 0.0, Test_Poll, NULL };
static void Test_Poll(void *unused)
{
	(void)unused; // using the unused for compiler because it asked nicely
	struct qsockaddr clientaddr;
	s32 control;
	s32 len;
	s8 name[32];
	s8 address[64];
	s32 colors;
	s32 frags;
	s32 connectTime;
	net_landriverlevel = testDriver;
	while (1) {
		len =
		    dfunc.Read(testSocket, net_message.data,
			       net_message.maxsize, &clientaddr);
		if (len < (s32)sizeof(s32))
			break;
		net_message.cursize = len;
		MSG_BeginReading();
		control = BigLong(*((s32 *)net_message.data));
		MSG_ReadLong();
		if (control == -1)
			break;
		if ((control & (~NETFLAG_LENGTH_MASK)) != (s32)NETFLAG_CTL)
			break;
		if ((control & NETFLAG_LENGTH_MASK) != len)
			break;
		if (MSG_ReadByte() != CCREP_PLAYER_INFO)
			Sys_Error
			    ("Unexpected response to Player Info request\n");
		MSG_ReadByte();	/* playerNumber */
		Q_strcpy(name, MSG_ReadString());
		colors = MSG_ReadLong();
		frags = MSG_ReadLong();
		connectTime = MSG_ReadLong();
		Q_strcpy(address, MSG_ReadString());
		Con_Printf("%s\n  frags:%3i  colors:%d %d  time:%d\n  %s\n",
			   name, frags, colors >> 4, colors & 0x0f,
			   connectTime / 60, address);
	}
	testPollCount--;
	if (testPollCount) {
		SchedulePollProcedure(&testPollProcedure, 0.1);
	} else {
		dfunc.Close_Socket(testSocket);
		testInProgress = 0;
	}
}

static void Test_f()
{
	const s8 *host;
	s32 n;
	s32 maxusers = MAX_SCOREBOARD;
	struct qsockaddr sendaddr;
	if (testInProgress)
		return;
	host = Strip_Port(Cmd_Argv(1));
	if (host && hostCacheCount) {
		for (n = 0; n < hostCacheCount; n++) {
			if (q_strcasecmp(host, hostcache[n].name) == 0) {
				if (hostcache[n].driver != myDriverLevel)
					continue;
				net_landriverlevel = hostcache[n].ldriver;
				maxusers = hostcache[n].maxusers;
				Q_memcpy(&sendaddr, &hostcache[n].addr,
					 sizeof(struct qsockaddr));
				break;
			}
		}
		if (n < hostCacheCount)
			goto JustDoIt;
	}
	for (net_landriverlevel = 0; net_landriverlevel < net_numlandrivers;
	     net_landriverlevel++) {
		if (!net_landrivers[net_landriverlevel].initialized)
			continue;
		// see if we can resolve the host name
		if (dfunc.GetAddrFromName(host, &sendaddr) != -1)
			break;
	}
	if (net_landriverlevel == net_numlandrivers) {
		Con_Printf("Could not resolve %s\n", host);
		return;
	}
JustDoIt:
	testSocket = dfunc.Open_Socket(0);
	if (testSocket == INVALID_SOCKET)
		return;
	testInProgress = 1;
	testPollCount = 20;
	testDriver = net_landriverlevel;
	for (n = 0; n < maxusers; n++) {
		SZ_Clear(&net_message);
		// save space for the header, filled in later
		MSG_WriteLong(&net_message, 0);
		MSG_WriteByte(&net_message, CCREQ_PLAYER_INFO);
		MSG_WriteByte(&net_message, n);
		*((s32 *)net_message.data) =
		    BigLong(NETFLAG_CTL |
			    (net_message.cursize & NETFLAG_LENGTH_MASK));
		dfunc.Write(testSocket, net_message.data, net_message.cursize,
			    &sendaddr);
	}
	SZ_Clear(&net_message);
	SchedulePollProcedure(&testPollProcedure, 0.1);
}

static bool test2InProgress = 0;
static s32 test2Driver;
static sys_socket_t test2Socket;
static void Test2_Poll(void *);
static PollProcedure test2PollProcedure = { NULL, 0.0, Test2_Poll, NULL };
static void Test2_Poll(SDL_UNUSED void *unused)
{
	struct qsockaddr clientaddr;
	s32 control;
	s32 len;
	s8 name[256];
	s8 value[256];
	net_landriverlevel = test2Driver;
	name[0] = 0;
	len =
	    dfunc.Read(test2Socket, net_message.data, net_message.maxsize,
		       &clientaddr);
	if (len < (s32)sizeof(s32))
		goto Reschedule;
	net_message.cursize = len;
	MSG_BeginReading();
	control = BigLong(*((s32 *)net_message.data));
	MSG_ReadLong();
	if (control == -1)
		goto Error;
	if ((control & (~NETFLAG_LENGTH_MASK)) != (s32)NETFLAG_CTL)
		goto Error;
	if ((control & NETFLAG_LENGTH_MASK) != len)
		goto Error;
	if (MSG_ReadByte() != CCREP_RULE_INFO)
		goto Error;
	Q_strcpy(name, MSG_ReadString());
	if (name[0] == 0)
		goto Done;
	Q_strcpy(value, MSG_ReadString());
	Con_Printf("%-16.16s  %-16.16s\n", name, value);
	SZ_Clear(&net_message);
	// save space for the header, filled in later
	MSG_WriteLong(&net_message, 0);
	MSG_WriteByte(&net_message, CCREQ_RULE_INFO);
	MSG_WriteString(&net_message, name);
	*((s32 *)net_message.data) =
	    BigLong(NETFLAG_CTL | (net_message.cursize & NETFLAG_LENGTH_MASK));
	dfunc.Write(test2Socket, net_message.data, net_message.cursize,
		    &clientaddr);
	SZ_Clear(&net_message);
Reschedule:
	SchedulePollProcedure(&test2PollProcedure, 0.05);
	return;
Error:
	Con_Printf("Unexpected response to Rule Info request\n");
Done:
	dfunc.Close_Socket(test2Socket);
	test2InProgress = 0;
	return;
}

static void Test2_f()
{
	const s8 *host;
	s32 n;
	struct qsockaddr sendaddr;
	if (test2InProgress)
		return;
	host = Strip_Port(Cmd_Argv(1));
	if (host && hostCacheCount) {
		for (n = 0; n < hostCacheCount; n++) {
			if (q_strcasecmp(host, hostcache[n].name) == 0) {
				if (hostcache[n].driver != myDriverLevel)
					continue;
				net_landriverlevel = hostcache[n].ldriver;
				Q_memcpy(&sendaddr, &hostcache[n].addr,
					 sizeof(struct qsockaddr));
				break;
			}
		}
		if (n < hostCacheCount)
			goto JustDoIt;
	}
	for (net_landriverlevel = 0; net_landriverlevel < net_numlandrivers;
	     net_landriverlevel++) {
		if (!net_landrivers[net_landriverlevel].initialized)
			continue;
		// see if we can resolve the host name
		if (dfunc.GetAddrFromName(host, &sendaddr) != -1)
			break;
	}
	if (net_landriverlevel == net_numlandrivers) {
		Con_Printf("Could not resolve %s\n", host);
		return;
	}
JustDoIt:
	test2Socket = dfunc.Open_Socket(0);
	if (test2Socket == INVALID_SOCKET)
		return;
	test2InProgress = 1;
	test2Driver = net_landriverlevel;
	SZ_Clear(&net_message);
	// save space for the header, filled in later
	MSG_WriteLong(&net_message, 0);
	MSG_WriteByte(&net_message, CCREQ_RULE_INFO);
	MSG_WriteString(&net_message, "");
	*((s32 *)net_message.data) =
	    BigLong(NETFLAG_CTL | (net_message.cursize & NETFLAG_LENGTH_MASK));
	dfunc.Write(test2Socket, net_message.data, net_message.cursize,
		    &sendaddr);
	SZ_Clear(&net_message);
	SchedulePollProcedure(&test2PollProcedure, 0.05);
}

s32 Datagram_Init()
{
	s32 i, num_inited;
	sys_socket_t csock;
#ifdef BAN_TEST
	banAddr.s_addr = INADDR_ANY;
	banMask.s_addr = INADDR_NONE;
#endif
	myDriverLevel = net_driverlevel;
	Cmd_AddCommand("net_stats", NET_Stats_f);
	s32 safemode = 0;
	if (safemode || COM_CheckParm("-nolan"))
		return -1;
	num_inited = 0;
	for (i = 0; i < net_numlandrivers; i++) {
		csock = net_landrivers[i].Init();
		if (csock == INVALID_SOCKET)
			continue;
		net_landrivers[i].initialized = 1;
		net_landrivers[i].controlSock = csock;
		num_inited++;
	}
	if (num_inited == 0)
		return -1;
#ifdef BAN_TEST
	Cmd_AddCommand("ban", NET_Ban_f);
#endif
	Cmd_AddCommand("test", Test_f);
	Cmd_AddCommand("test2", Test2_f);
	return 0;
}

void Datagram_Shutdown()
{
	s32 i;
// shutdown the lan drivers
	for (i = 0; i < net_numlandrivers; i++) {
		if (net_landrivers[i].initialized) {
			net_landrivers[i].Shutdown();
			net_landrivers[i].initialized = 0;
		}
	}
}

void Datagram_Close(qsocket_t *sock)
{
	sfunc.Close_Socket(sock->socket);
}

void Datagram_Listen(bool state)
{
	s32 i;
	for (i = 0; i < net_numlandrivers; i++) {
		if (net_landrivers[i].initialized)
			net_landrivers[i].Listen(state);
	}
}

static qsocket_t *_Datagram_CheckNewConnections()
{
	struct qsockaddr clientaddr;
	struct qsockaddr newaddr;
	sys_socket_t newsock;
	sys_socket_t acceptsock;
	qsocket_t *sock;
	qsocket_t *s;
	s32 len;
	s32 command;
	s32 control;
	s32 ret;
	acceptsock = dfunc.CheckNewConnections();
	if (acceptsock == INVALID_SOCKET)
		return NULL;
	SZ_Clear(&net_message);
	len =
	    dfunc.Read(acceptsock, net_message.data, net_message.maxsize,
		       &clientaddr);
	if (len < (s32)sizeof(s32))
		return NULL;
	net_message.cursize = len;
	MSG_BeginReading();
	control = BigLong(*((s32 *)net_message.data));
	MSG_ReadLong();
	if (control == -1)
		return NULL;
	if ((control & (~NETFLAG_LENGTH_MASK)) != (s32)NETFLAG_CTL)
		return NULL;
	if ((control & NETFLAG_LENGTH_MASK) != len)
		return NULL;
	command = MSG_ReadByte();
	if (command == CCREQ_SERVER_INFO) {
		if (Q_strcmp(MSG_ReadString(), "QUAKE") != 0)
			return NULL;
		SZ_Clear(&net_message);
		// save space for the header, filled in later
		MSG_WriteLong(&net_message, 0);
		MSG_WriteByte(&net_message, CCREP_SERVER_INFO);
		dfunc.GetSocketAddr(acceptsock, &newaddr);
		MSG_WriteString(&net_message, (s8*)dfunc.AddrToString(&newaddr));
		MSG_WriteString(&net_message, hostname.string);
		MSG_WriteString(&net_message, sv.name);
		MSG_WriteByte(&net_message, net_activeconnections);
		MSG_WriteByte(&net_message, svs.maxclients);
		MSG_WriteByte(&net_message, NET_PROTOCOL_VERSION);
		*((s32 *)net_message.data) =
		    BigLong(NETFLAG_CTL |
			    (net_message.cursize & NETFLAG_LENGTH_MASK));
		dfunc.Write(acceptsock, net_message.data, net_message.cursize,
			    &clientaddr);
		SZ_Clear(&net_message);
		return NULL;
	}
	if (command == CCREQ_PLAYER_INFO) {
		s32 playerNumber;
		s32 activeNumber;
		s32 clientNumber;
		client_t *client;
		playerNumber = MSG_ReadByte();
		activeNumber = -1;
		for (clientNumber = 0, client = svs.clients;
		     clientNumber < svs.maxclients; clientNumber++, client++) {
			if (client->active) {
				activeNumber++;
				if (activeNumber == playerNumber)
					break;
			}
		}
		if (clientNumber == svs.maxclients)
			return NULL;
		SZ_Clear(&net_message);
		// save space for the header, filled in later
		MSG_WriteLong(&net_message, 0);
		MSG_WriteByte(&net_message, CCREP_PLAYER_INFO);
		MSG_WriteByte(&net_message, playerNumber);
		MSG_WriteString(&net_message, client->name);
		MSG_WriteLong(&net_message, client->colors);
		MSG_WriteLong(&net_message, (s32)client->edict->v.frags);
		MSG_WriteLong(&net_message,
			      (s32)(net_time -
				    client->netconnection->connecttime));
		MSG_WriteString(&net_message, client->netconnection->address);
		*((s32 *)net_message.data) =
		    BigLong(NETFLAG_CTL |
			    (net_message.cursize & NETFLAG_LENGTH_MASK));
		dfunc.Write(acceptsock, net_message.data, net_message.cursize,
			    &clientaddr);
		SZ_Clear(&net_message);
		return NULL;
	}
	if (command == CCREQ_RULE_INFO) {
		const s8 *prevCvarName;
		cvar_t *var;
		// find the search start location
		prevCvarName = MSG_ReadString();
		var = Cvar_FindVar(prevCvarName);
		// send the response
		SZ_Clear(&net_message);
		// save space for the header, filled in later
		MSG_WriteLong(&net_message, 0);
		MSG_WriteByte(&net_message, CCREP_RULE_INFO);
		if (var) {
			MSG_WriteString(&net_message, var->name);
			MSG_WriteString(&net_message, var->string);
		}
		*((s32 *)net_message.data) =
		    BigLong(NETFLAG_CTL |
			    (net_message.cursize & NETFLAG_LENGTH_MASK));
		dfunc.Write(acceptsock, net_message.data, net_message.cursize,
			    &clientaddr);
		SZ_Clear(&net_message);
		return NULL;
	}
	if (command != CCREQ_CONNECT)
		return NULL;
	if (Q_strcmp(MSG_ReadString(), "QUAKE") != 0)
		return NULL;
	if (MSG_ReadByte() != NET_PROTOCOL_VERSION) {
		SZ_Clear(&net_message);
		// save space for the header, filled in later
		MSG_WriteLong(&net_message, 0);
		MSG_WriteByte(&net_message, CCREP_REJECT);
		MSG_WriteString(&net_message, "Incompatible version.\n");
		*((s32 *)net_message.data) =
		    BigLong(NETFLAG_CTL |
			    (net_message.cursize & NETFLAG_LENGTH_MASK));
		dfunc.Write(acceptsock, net_message.data, net_message.cursize,
			    &clientaddr);
		SZ_Clear(&net_message);
		return NULL;
	}
#ifdef BAN_TEST
	// check for a ban
	if (clientaddr.qsa_family == AF_INET) {
		in_addr_t testAddr;
		testAddr = ((struct sockaddr_in *)&clientaddr)->sin_addr.s_addr;
		if ((testAddr & banMask.s_addr) == banAddr.s_addr) {
			SZ_Clear(&net_message);
			// save space for the header, filled in later
			MSG_WriteLong(&net_message, 0);
			MSG_WriteByte(&net_message, CCREP_REJECT);
			MSG_WriteString(&net_message,
					"You have been banned.\n");
			*((s32 *)net_message.data) =
			    BigLong(NETFLAG_CTL |
				    (net_message.
				     cursize & NETFLAG_LENGTH_MASK));
			dfunc.Write(acceptsock, net_message.data,
				    net_message.cursize, &clientaddr);
			SZ_Clear(&net_message);
			return NULL;
		}
	}
#endif
	// see if this guy is already connected
	for (s = net_activeSockets; s; s = s->next) {
		if (s->driver != net_driverlevel)
			continue;
		ret = dfunc.AddrCompare(&clientaddr, &s->addr);
		if (ret >= 0) {
			// is this a duplicate connection reqeust?
			if (ret == 0 && net_time - s->connecttime < 2.0) {
				// yes, so send a duplicate reply
				SZ_Clear(&net_message);
				// save space for the header, filled in later
				MSG_WriteLong(&net_message, 0);
				MSG_WriteByte(&net_message, CCREP_ACCEPT);
				dfunc.GetSocketAddr(s->socket, &newaddr);
				MSG_WriteLong(&net_message,
					      dfunc.GetSocketPort(&newaddr));
				*((s32 *)net_message.data) =
				    BigLong(NETFLAG_CTL |
					    (net_message.
					     cursize & NETFLAG_LENGTH_MASK));
				dfunc.Write(acceptsock, net_message.data,
					    net_message.cursize, &clientaddr);
				SZ_Clear(&net_message);
				return NULL;
			}
			// it's somebody coming back in from a crash/disconnect
			// so close the old qsocket and let their retry get them back in
			NET_Close(s);
			return NULL;
		}
	}
	// allocate a QSocket
	sock = NET_NewQSocket();
	if (sock == NULL)	// no room; try to let him know
	{
		SZ_Clear(&net_message);
		// save space for the header, filled in later
		MSG_WriteLong(&net_message, 0);
		MSG_WriteByte(&net_message, CCREP_REJECT);
		MSG_WriteString(&net_message, "Server is full.\n");
		*((s32 *)net_message.data) =
		    BigLong(NETFLAG_CTL |
			    (net_message.cursize & NETFLAG_LENGTH_MASK));
		dfunc.Write(acceptsock, net_message.data, net_message.cursize,
			    &clientaddr);
		SZ_Clear(&net_message);
		return NULL;
	}
	// allocate a network socket
	newsock = dfunc.Open_Socket(0);
	if (newsock == INVALID_SOCKET) {
		NET_FreeQSocket(sock);
		return NULL;
	}
	// connect to the client
	if (dfunc.Connect(newsock, &clientaddr) == -1) {
		dfunc.Close_Socket(newsock);
		NET_FreeQSocket(sock);
		return NULL;
	}
	// everything is allocated, just fill in the details
	sock->socket = newsock;
	sock->landriver = net_landriverlevel;
	sock->addr = clientaddr;
	Q_strcpy(sock->address, dfunc.AddrToString(&clientaddr));
	// send him back the info about the server connection he has been allocated
	SZ_Clear(&net_message);
	// save space for the header, filled in later
	MSG_WriteLong(&net_message, 0);
	MSG_WriteByte(&net_message, CCREP_ACCEPT);
	dfunc.GetSocketAddr(newsock, &newaddr);
	MSG_WriteLong(&net_message, dfunc.GetSocketPort(&newaddr));
	*((s32 *)net_message.data) =
	    BigLong(NETFLAG_CTL | (net_message.cursize & NETFLAG_LENGTH_MASK));
	dfunc.Write(acceptsock, net_message.data, net_message.cursize,
		    &clientaddr);
	SZ_Clear(&net_message);
	return sock;
}

qsocket_t *Datagram_CheckNewConnections()
{
	qsocket_t *ret = NULL;
	for (net_landriverlevel = 0; net_landriverlevel < net_numlandrivers;
	     net_landriverlevel++) {
		if (net_landrivers[net_landriverlevel].initialized) {
			if ((ret = _Datagram_CheckNewConnections()) != NULL)
				break;
		}
	}
	return ret;
}

static void _Datagram_SearchForHosts(bool xmit)
{
	s32 ret;
	s32 n;
	s32 i;
	struct qsockaddr readaddr;
	struct qsockaddr myaddr;
	s32 control;
	dfunc.GetSocketAddr(dfunc.controlSock, &myaddr);
	if (xmit) {
		SZ_Clear(&net_message);
		// save space for the header, filled in later
		MSG_WriteLong(&net_message, 0);
		MSG_WriteByte(&net_message, CCREQ_SERVER_INFO);
		MSG_WriteString(&net_message, "QUAKE");
		MSG_WriteByte(&net_message, NET_PROTOCOL_VERSION);
		*((s32 *)net_message.data) =
		    BigLong(NETFLAG_CTL |
			    (net_message.cursize & NETFLAG_LENGTH_MASK));
		dfunc.Broadcast(dfunc.controlSock, net_message.data,
				net_message.cursize);
		SZ_Clear(&net_message);
	}
	while ((ret =
		dfunc.Read(dfunc.controlSock, net_message.data,
			   net_message.maxsize, &readaddr)) > 0) {
		if (ret < (s32)sizeof(s32))
			continue;
		net_message.cursize = ret;
		// don't answer our own query
		if (dfunc.AddrCompare(&readaddr, &myaddr) >= 0)
			continue;
		// is the cache full?
		if (hostCacheCount == HOSTCACHESIZE)
			continue;
		MSG_BeginReading();
		control = BigLong(*((s32 *)net_message.data));
		MSG_ReadLong();
		if (control == -1)
			continue;
		if ((control & (~NETFLAG_LENGTH_MASK)) != (s32)NETFLAG_CTL)
			continue;
		if ((control & NETFLAG_LENGTH_MASK) != ret)
			continue;
		if (MSG_ReadByte() != CCREP_SERVER_INFO)
			continue;
		dfunc.GetAddrFromName(MSG_ReadString(), &readaddr);
		// search the cache for this server
		for (n = 0; n < hostCacheCount; n++) {
			if (dfunc.AddrCompare(&readaddr, &hostcache[n].addr) ==
			    0)
				break;
		}
		// is it already there?
		if (n < hostCacheCount)
			continue;
		// add it
		hostCacheCount++;
		Q_strcpy(hostcache[n].name, MSG_ReadString());
		Q_strcpy(hostcache[n].map, MSG_ReadString());
		hostcache[n].users = MSG_ReadByte();
		hostcache[n].maxusers = MSG_ReadByte();
		if (MSG_ReadByte() != NET_PROTOCOL_VERSION) {
			Q_strcpy(hostcache[n].cname, hostcache[n].name);
			hostcache[n].cname[14] = 0;
			Q_strcpy(hostcache[n].name, "*");
			Q_strcat(hostcache[n].name, hostcache[n].cname);
		}
		Q_memcpy(&hostcache[n].addr, &readaddr,
			 sizeof(struct qsockaddr));
		hostcache[n].driver = net_driverlevel;
		hostcache[n].ldriver = net_landriverlevel;
		Q_strcpy(hostcache[n].cname, dfunc.AddrToString(&readaddr));
		// check for a name conflict
		for (i = 0; i < hostCacheCount; i++) {
			if (i == n)
				continue;
			if (q_strcasecmp(hostcache[n].name, hostcache[i].name)
			    == 0) {
				i = Q_strlen(hostcache[n].name);
				if (i < 15 && hostcache[n].name[i - 1] > '8') {
					hostcache[n].name[i] = '0';
					hostcache[n].name[i + 1] = 0;
				} else
					hostcache[n].name[i - 1]++;
				i = -1;
			}
		}
	}
}

void Datagram_SearchForHosts(bool xmit)
{
	for (net_landriverlevel = 0; net_landriverlevel < net_numlandrivers;
	     net_landriverlevel++) {
		if (hostCacheCount == HOSTCACHESIZE)
			break;
		if (net_landrivers[net_landriverlevel].initialized)
			_Datagram_SearchForHosts(xmit);
	}
}

static qsocket_t *_Datagram_Connect(const s8 *host)
{
	struct qsockaddr sendaddr;
	struct qsockaddr readaddr;
	qsocket_t *sock;
	sys_socket_t newsock;
	s32 ret;
	s32 reps;
	f64 start_time;
	s32 control;
	const s8 *reason;
	// see if we can resolve the host name
	if (dfunc.GetAddrFromName(host, &sendaddr) == -1) {
		Con_Printf("Could not resolve %s\n", host);
		return NULL;
	}
	newsock = dfunc.Open_Socket(0);
	if (newsock == INVALID_SOCKET)
		return NULL;
	sock = NET_NewQSocket();
	if (sock == NULL)
		goto ErrorReturn2;
	sock->socket = newsock;
	sock->landriver = net_landriverlevel;
	// connect to the host
	if (dfunc.Connect(newsock, &sendaddr) == -1)
		goto ErrorReturn;
	// send the connection request
	Con_Printf("trying...\n");
	SCR_UpdateScreen();
	start_time = net_time;
	for (reps = 0; reps < 3; reps++) {
		SZ_Clear(&net_message);
		// save space for the header, filled in later
		MSG_WriteLong(&net_message, 0);
		MSG_WriteByte(&net_message, CCREQ_CONNECT);
		MSG_WriteString(&net_message, "QUAKE");
		MSG_WriteByte(&net_message, NET_PROTOCOL_VERSION);
		*((s32 *)net_message.data) =
		    BigLong(NETFLAG_CTL |
			    (net_message.cursize & NETFLAG_LENGTH_MASK));
		dfunc.Write(newsock, net_message.data, net_message.cursize,
			    &sendaddr);
		SZ_Clear(&net_message);
		do {
			ret =
			    dfunc.Read(newsock, net_message.data,
				       net_message.maxsize, &readaddr);
			// if we got something, validate it
			if (ret > 0) {
				// is it from the right place?
				if (sfunc.AddrCompare(&readaddr, &sendaddr) !=
				    0) {
					Con_Printf("wrong reply address\n");
					Con_Printf("Expected: %s | %s\n",
						   dfunc.
						   AddrToString(&sendaddr),
						   StrAddr(&sendaddr));
					Con_Printf("Received: %s | %s\n",
						   dfunc.
						   AddrToString(&readaddr),
						   StrAddr(&readaddr));
					SCR_UpdateScreen();
					ret = 0;
					continue;
				}
				if (ret < (s32)sizeof(s32)) {
					ret = 0;
					continue;
				}
				net_message.cursize = ret;
				MSG_BeginReading();
				control = BigLong(*((s32 *)net_message.data));
				MSG_ReadLong();
				if (control == -1) {
					ret = 0;
					continue;
				}
				if ((control & (~NETFLAG_LENGTH_MASK)) !=
				    (s32)NETFLAG_CTL) {
					ret = 0;
					continue;
				}
				if ((control & NETFLAG_LENGTH_MASK) != ret) {
					ret = 0;
					continue;
				}
			}
		}
		while (ret == 0 && (SetNetTime() - start_time) < 2.5);
		if (ret)
			break;
		Con_Printf("still trying...\n");
		SCR_UpdateScreen();
		start_time = SetNetTime();
	}
	if (ret == 0) {
		reason = "No Response";
		Con_Printf("%s\n", reason);
		Q_strcpy(m_return_reason, reason);
		goto ErrorReturn;
	}
	if (ret == -1) {
		reason = "Network Error";
		Con_Printf("%s\n", reason);
		Q_strcpy(m_return_reason, reason);
		goto ErrorReturn;
	}
	ret = MSG_ReadByte();
	if (ret == CCREP_REJECT) {
		reason = MSG_ReadString();
		Con_Printf("%s\n", reason);
		Q_strncpy(m_return_reason, reason, sizeof(m_return_reason));
		goto ErrorReturn;
	}
	if (ret == CCREP_ACCEPT) {
		Q_memcpy(&sock->addr, &sendaddr, sizeof(struct qsockaddr));
		dfunc.SetSocketPort(&sock->addr, MSG_ReadLong());
	} else {
		reason = "Bad Response";
		Con_Printf("%s\n", reason);
		Q_strcpy(m_return_reason, reason);
		goto ErrorReturn;
	}
	dfunc.GetNameFromAddr(&sendaddr, sock->address);
	Con_Printf("Connection accepted\n");
	sock->lastMessageTime = SetNetTime();
	// switch the connection to the specified address
	if (dfunc.Connect(newsock, &sock->addr) == -1) {
		reason = "Connect to Game failed";
		Con_Printf("%s\n", reason);
		Q_strcpy(m_return_reason, reason);
		goto ErrorReturn;
	}
	m_return_onerror = 0;
	return sock;
ErrorReturn:
	NET_FreeQSocket(sock);
ErrorReturn2:
	dfunc.Close_Socket(newsock);
	if (m_return_onerror) {
		//IN_DeactivateForMenu();
		key_dest = key_menu;
		//m_state = m_return_state;
		m_return_onerror = 0;
	}
	return NULL;
}

qsocket_t *Datagram_Connect(const s8 *host)
{
	qsocket_t *ret = NULL;
	host = Strip_Port(host);
	for (net_landriverlevel = 0; net_landriverlevel < net_numlandrivers;
	     net_landriverlevel++) {
		if (net_landrivers[net_landriverlevel].initialized) {
			if ((ret = _Datagram_Connect(host)) != NULL)
				break;
		}
	}
	return ret;
}
#undef sfunc
#undef dfunc
