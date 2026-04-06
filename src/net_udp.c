// Copyright (C) 1996-2001 Id Software, Inc.
// Copyright (C) 2007-2008 Kristian Duske
// Copyright (C) 2010-2014 QuakeSpasm developers
// GPLv3 See LICENSE for details.

#include "quakedef.h"

static sys_socket_t net_acceptsocket = INVALID_SOCKET;	// socket for fielding new connections
static sys_socket_t net_controlsocket;
static sys_socket_t net_broadcastsocket = 0;
static struct sockaddr_in broadcastaddr;
static in_addr_t myAddr;

sys_socket_t UDP_Init()
{
	s8 buff[MAXHOSTNAMELEN];
	struct hostent *local;
	struct qsockaddr addr;
	if (COM_CheckParm("-noudp"))
		return INVALID_SOCKET;
	// determine my name & address
	myAddr = htonl(INADDR_LOOPBACK);
	if (gethostname(buff, MAXHOSTNAMELEN) != 0) {
		s32 err = SOCKETERRNO;
		Con_SafePrintf("UDP_Init: WARNING: gethostname failed (%s)\n",
			       socketerror(err));
	} else {
		buff[MAXHOSTNAMELEN - 1] = 0;
		if (!(local = gethostbyname(buff))) {
#ifdef _WIN32
			const s8 herrmsg[1024];
			FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER | 
				FORMAT_MESSAGE_FROM_SYSTEM | 
				FORMAT_MESSAGE_IGNORE_INSERTS, 
				NULL, 
				h_errno, 
				MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), 
				(LPWSTR)&herrmsg, 0, NULL);
			Con_SafePrintf
			("UDP_Init: WARNING: gethostbyname failed (%s)\n",
				herrmsg);

#else
			Con_SafePrintf
			    ("UDP_Init: WARNING: gethostbyname failed (%s)\n",
			     hstrerror(h_errno));
#endif
		} else if (local->h_addrtype != AF_INET) {
			Con_SafePrintf
			    ("UDP_Init: address from gethostbyname not IPv4\n");
		} else {
			s32 i = COM_CheckParm("-ip");
			if (i) {
				if (i < com_argc - 1) {
					myAddr = inet_addr(com_argv[i + 1]);
					if (myAddr == INADDR_NONE)
						Sys_Error ("%s is not a valid IP address", com_argv[i + 1]);
					strcpy(my_tcpip_address, com_argv[i+1]);
				} else {
					Sys_Error ("NET_Init: you must specify an IP address after -ip");
				}
			} else {
				myAddr = *(in_addr_t *) local->h_addr_list[0];
			}
		}
	}
	if ((net_controlsocket = UDP_OpenSocket(0)) == INVALID_SOCKET) {
		Con_SafePrintf
		    ("UDP_Init: Unable to open control socket, UDP disabled\n");
		return INVALID_SOCKET;
	}
	broadcastaddr.sin_family = AF_INET;
	broadcastaddr.sin_addr.s_addr = INADDR_BROADCAST;
	broadcastaddr.sin_port = htons((u16)net_hostport);
	UDP_GetSocketAddr(net_controlsocket, &addr);
	strcpy(my_tcpip_address, UDP_AddrToString(&addr));
	s8 *tst = strrchr(my_tcpip_address, ':');
	if (tst)
		*tst = 0;
	Con_SafePrintf("UDP Initialized\n");
	tcpipAvailable = 1;
	return net_controlsocket;
}

void UDP_Shutdown()
{
	UDP_Listen(0);
	UDP_CloseSocket(net_controlsocket);
}

void UDP_Listen(bool state)
{
	// enable listening
	if (state) {
		if (net_acceptsocket != INVALID_SOCKET)
			return;
		if ((net_acceptsocket =
		     UDP_OpenSocket(net_hostport)) == INVALID_SOCKET)
			Sys_Error("UDP_Listen: Unable to open accept socket");
		return;
	}
	// disable listening
	if (net_acceptsocket == INVALID_SOCKET)
		return;
	UDP_CloseSocket(net_acceptsocket);
	net_acceptsocket = INVALID_SOCKET;
}

sys_socket_t UDP_OpenSocket(s32 port)
{
	sys_socket_t newsocket;
	struct sockaddr_in address;
	s32 _true = 1;
	s32 err;
	if ((newsocket =
	     socket(PF_INET, SOCK_DGRAM, IPPROTO_UDP)) == INVALID_SOCKET) {
		err = SOCKETERRNO;
		Con_SafePrintf("UDP_OpenSocket: %s\n", socketerror(err));
		return INVALID_SOCKET;
	}
	if (ioctlsocket(newsocket, FIONBIO, &_true) == SOCKET_ERROR)
		goto ErrorReturn;
	memset(&address, 0, sizeof(struct sockaddr_in));
	address.sin_family = AF_INET;
	address.sin_addr.s_addr = INADDR_ANY;
	address.sin_port = htons((u16)port);
	if (bind(newsocket, (struct sockaddr *)&address, sizeof(address)) == 0)
		return newsocket;
ErrorReturn:
	err = SOCKETERRNO;
	Con_SafePrintf("UDP_OpenSocket: %s\n", socketerror(err));
	UDP_CloseSocket(newsocket);
	return INVALID_SOCKET;
}

s32 UDP_CloseSocket(sys_socket_t socketid)
{
	if (socketid == net_broadcastsocket)
		net_broadcastsocket = 0;
	return closesocket(socketid);
}
/*
PartialIPAddress
this lets you type only as much of the net address as required, using
the local network components to fill in the rest
*/
static s32 PartialIPAddress(const s8 *in, struct qsockaddr *hostaddr)
{
	s8 buff[256];
	buff[0] = '.';
	s8 *b = buff;
	strcpy(buff + 1, in);
	if (buff[1] == '.')
		b++;
	s32 addr = 0;
	s32 mask = -1;
	while (*b == '.') {
		b++;
		s32 num = 0;
		s32 run = 0;
		while (!(*b < '0' || *b > '9')) {
			num = num * 10 + *b++ - '0';
			if (++run > 3)
				return -1;
		}
		if ((*b < '0' || *b > '9') && *b != '.' && *b != ':' && *b != 0)
			return -1;
		if (num < 0 || num > 255)
			return -1;
		mask <<= 8;
		addr = (addr << 8) + num;
	}
	s32 port = net_hostport;
	if (*b++ == ':')
		port = atoi(b);
	hostaddr->qsa_family = AF_INET;
	((struct sockaddr_in*)hostaddr)->sin_port = htons((u16)port);
	((struct sockaddr_in*)hostaddr)->sin_addr.s_addr =
	    (myAddr & htonl(mask)) | htonl(addr);
	return 0;
}

s32 UDP_Connect(sys_socket_t socketid, struct qsockaddr *addr)
{
	(void)socketid;
	(void)addr;
	return 0;
}

sys_socket_t UDP_CheckNewConnections()
{
	s32 available;
	struct sockaddr_in from;
	socklen_t fromlen;
	s8 buff[1];
	if (net_acceptsocket == INVALID_SOCKET)
		return INVALID_SOCKET;
	if (ioctlsocket(net_acceptsocket, FIONREAD, &available) == -1) {
		s32 err = SOCKETERRNO;
		Sys_Error("UDP: ioctlsocket (FIONREAD) failed (%s)",
			  socketerror(err));
	}
	if (available)
		return net_acceptsocket;
	// quietly absorb empty packets
	recvfrom(net_acceptsocket, buff, 0, 0,(struct sockaddr*)&from,&fromlen);
	return INVALID_SOCKET;
}

s32 UDP_Read(sys_socket_t socketid, u8 *buf, s32 len, struct qsockaddr *addr)
{
	socklen_t addrlen = sizeof(struct qsockaddr);
	s32 ret = recvfrom(socketid,buf,len,0,(struct sockaddr*)addr,&addrlen);
	if (ret == SOCKET_ERROR) {
		s32 err = SOCKETERRNO;
		if (err == NET_EWOULDBLOCK || err == NET_ECONNREFUSED)
			return 0;
		Con_SafePrintf("UDP_Read, recvfrom: %s\n", socketerror(err));
	}
	return ret;
}

static s32 UDP_MakeSocketBroadcastCapable(sys_socket_t socketid)
{
	s32 i = 1;
	// make this socket broadcast capable
	if (setsockopt(socketid, SOL_SOCKET, SO_BROADCAST, (s8*)&i, sizeof(i))
	    == SOCKET_ERROR) {
		s32 err = SOCKETERRNO;
		Con_SafePrintf("UDP, setsockopt: %s\n", socketerror(err));
		return -1;
	}
	net_broadcastsocket = socketid;
	return 0;
}

s32 UDP_Broadcast(sys_socket_t socketid, u8 *buf, s32 len)
{
	if (socketid != net_broadcastsocket) {
		if (net_broadcastsocket != 0)
			Sys_Error
			    ("Attempted to use multiple broadcasts sockets");
		s32 ret = UDP_MakeSocketBroadcastCapable(socketid);
		if (ret == -1) {
			Con_Printf("Unable to make socket broadcast capable\n");
			return ret;
		}
	}
	return UDP_Write(socketid, buf, len,
			 (struct qsockaddr *)&broadcastaddr);
}

s32 UDP_Write(sys_socket_t socketid, u8 *buf, s32 len, struct qsockaddr *addr)
{
	s32 ret = sendto(socketid, buf, len, 0, (struct sockaddr *)addr,
		     sizeof(struct qsockaddr));
	if (ret == SOCKET_ERROR) {
		s32 err = SOCKETERRNO;
		if (err == NET_EWOULDBLOCK)
			return 0;
		Con_SafePrintf("UDP_Write, sendto: %s\n", socketerror(err));
	}
	return ret;
}

const s8 *UDP_AddrToString(struct qsockaddr *addr)
{
	static s8 buffer[22];
	s32 haddr = ntohl(((struct sockaddr_in *)addr)->sin_addr.s_addr);
	snprintf(buffer, sizeof(buffer), "%d.%d.%d.%d:%d", (haddr >> 24) & 0xff,
		 (haddr >> 16) & 0xff, (haddr >> 8) & 0xff, haddr & 0xff,
		 ntohs(((struct sockaddr_in *)addr)->sin_port));
	return buffer;
}

s32 UDP_StringToAddr(const s8 *string, struct qsockaddr *addr)
{
	s32 ha1, ha2, ha3, ha4, hp, ipaddr;
	sscanf(string, "%d.%d.%d.%d:%d", &ha1, &ha2, &ha3, &ha4, &hp);
	ipaddr = (ha1 << 24) | (ha2 << 16) | (ha3 << 8) | ha4;
	addr->qsa_family = AF_INET;
	((struct sockaddr_in *)addr)->sin_addr.s_addr = htonl(ipaddr);
	((struct sockaddr_in *)addr)->sin_port = htons((u16)hp);
	return 0;
}

s32 UDP_GetSocketAddr(sys_socket_t socketid, struct qsockaddr *addr)
{
	socklen_t addrlen = sizeof(struct qsockaddr);
	in_addr_t a;
	memset(addr, 0, sizeof(struct qsockaddr));
	if (getsockname(socketid, (struct sockaddr *)addr, &addrlen) != 0)
		return -1;
	a = ((struct sockaddr_in *)addr)->sin_addr.s_addr;
	if (a == 0 || a == htonl(INADDR_LOOPBACK))
		((struct sockaddr_in *)addr)->sin_addr.s_addr = myAddr;
	return 0;
}

s32 UDP_GetNameFromAddr(struct qsockaddr *addr, s8 *name)
{
	struct hostent *hostentry = gethostbyaddr((s8*)&((struct sockaddr_in*)
			addr)->sin_addr, sizeof(struct in_addr), AF_INET);
	if (hostentry) {
		strncpy(name, (s8 *)hostentry->h_name, NET_NAMELEN - 1);
		return 0;
	}
	strcpy(name, UDP_AddrToString(addr));
	return 0;
}

s32 UDP_GetAddrFromName(const s8 *name, struct qsockaddr *addr)
{
	struct hostent *hostentry;
	if (name[0] >= '0' && name[0] <= '9')
		return PartialIPAddress(name, addr);
	hostentry = gethostbyname(name);
	if (!hostentry)
		return -1;
	addr->qsa_family = AF_INET;
	((struct sockaddr_in *)addr)->sin_port =
	    htons((u16)net_hostport);
	((struct sockaddr_in *)addr)->sin_addr.s_addr =
	    *(in_addr_t *) hostentry->h_addr_list[0];
	return 0;
}

s32 UDP_AddrCompare(struct qsockaddr *addr1, struct qsockaddr *addr2)
{
	if (addr1->qsa_family != addr2->qsa_family)
		return -1;
	if (((struct sockaddr_in *)addr1)->sin_addr.s_addr !=
	    ((struct sockaddr_in *)addr2)->sin_addr.s_addr)
		return -1;
	if (((struct sockaddr_in *)addr1)->sin_port !=
	    ((struct sockaddr_in *)addr2)->sin_port)
		return 1;
	return 0;
}

s32 UDP_GetSocketPort(struct qsockaddr *addr)
{
	return ntohs(((struct sockaddr_in *)addr)->sin_port);
}

s32 UDP_SetSocketPort(struct qsockaddr *addr, s32 port)
{
	((struct sockaddr_in *)addr)->sin_port = htons((u16)port);
	return 0;
}
