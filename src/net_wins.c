// Copyright (C) 1996-2001 Id Software, Inc.
// Copyright (C) 2010-2014 QuakeSpasm developers
// GPLv3 See LICENSE for details.

#include "quakedef.h"
#include "net_wins.h"

static sys_socket_t net_acceptsocket = INVALID_SOCKET;	// socket for fielding new connections
static sys_socket_t net_controlsocket;
static sys_socket_t net_broadcastsocket = 0;
static struct sockaddr_in broadcastaddr;
static in_addr_t myAddr;

s32 winsock_initialized = 0;
WSADATA winsockdata;

const s8 *__WSAE_StrError (s32 err)
{
	switch (err) {
	case 0:			return "No error";
	case WSAEINTR:		return "Interrupted system call";		/* 10004 */
	case WSAEBADF:		return "Bad file number";			/* 10009 */
	case WSAEACCES:		return "Permission denied";			/* 10013 */
	case WSAEFAULT:		return "Bad address";				/* 10014 */
	case WSAEINVAL:		return "Invalid argument (not bind)";		/* 10022 */
	case WSAEMFILE:		return "Too many open files";			/* 10024 */
	case WSAEWOULDBLOCK:	return "Operation would block";			/* 10035 */
	case WSAEINPROGRESS:	return "Operation now in progress";		/* 10036 */
	case WSAEALREADY:	return "Operation already in progress";		/* 10037 */
	case WSAENOTSOCK:	return "Socket operation on non-socket";	/* 10038 */
	case WSAEDESTADDRREQ:	return "Destination address required";		/* 10039 */
	case WSAEMSGSIZE:	return "Message too s64";			/* 10040 */
	case WSAEPROTOTYPE:	return "Protocol wrong type for socket";	/* 10041 */
	case WSAENOPROTOOPT:	return "Bad protocol option";			/* 10042 */
	case WSAEPROTONOSUPPORT: return "Protocol not supported";		/* 10043 */
	case WSAESOCKTNOSUPPORT: return "Socket type not supported";		/* 10044 */
	case WSAEOPNOTSUPP:	return "Operation not supported on socket";	/* 10045 */
	case WSAEPFNOSUPPORT:	return "Protocol family not supported";		/* 10046 */
	case WSAEAFNOSUPPORT:	return "Address family not supported by protocol family"; /* 10047 */
	case WSAEADDRINUSE:	return "Address already in use";		/* 10048 */
	case WSAEADDRNOTAVAIL:	return "Can't assign requested address";	/* 10049 */
	case WSAENETDOWN:	return "Network is down";			/* 10050 */
	case WSAENETUNREACH:	return "Network is unreachable";		/* 10051 */
	case WSAENETRESET:	return "Net dropped connection or reset";	/* 10052 */
	case WSAECONNABORTED:	return "Software caused connection abort";	/* 10053 */
	case WSAECONNRESET:	return "Connection reset by peer";		/* 10054 */
	case WSAENOBUFS:	return "No buffer space available";		/* 10055 */
	case WSAEISCONN:	return "Socket is already connected";		/* 10056 */
	case WSAENOTCONN:	return "Socket is not connected";		/* 10057 */
	case WSAESHUTDOWN:	return "Can't send after socket shutdown";	/* 10058 */
	case WSAETOOMANYREFS:	return "Too many references, can't splice";	/* 10059 */
	case WSAETIMEDOUT:	return "Connection timed out";			/* 10060 */
	case WSAECONNREFUSED:	return "Connection refused";			/* 10061 */
	case WSAELOOP:		return "Too many levels of symbolic links";	/* 10062 */
	case WSAENAMETOOLONG:	return "File name too s64";			/* 10063 */
	case WSAEHOSTDOWN:	return "Host is down";				/* 10064 */
	case WSAEHOSTUNREACH:	return "No Route to Host";			/* 10065 */
	case WSAENOTEMPTY:	return "Directory not empty";			/* 10066 */
	case WSAEPROCLIM:	return "Too many processes";			/* 10067 */
	case WSAEUSERS:		return "Too many users";			/* 10068 */
	case WSAEDQUOT:		return "Disc Quota Exceeded";			/* 10069 */
	case WSAESTALE:		return "Stale NFS file handle";			/* 10070 */
	case WSAEREMOTE:	return "Too many levels of remote in path";	/* 10071 */
	case WSAEDISCON:	return "Graceful shutdown in progress";		/* 10101 */
	case WSASYSNOTREADY:	return "Network SubSystem is unavailable";			/* 10091 */
	case WSAVERNOTSUPPORTED: return "WINSOCK DLL Version out of range";			/* 10092 */
	case WSANOTINITIALISED:	return "Successful WSASTARTUP not yet performed";		/* 10093 */
	case WSAHOST_NOT_FOUND:	return "Authoritative answer: Host not found";			/* 11001 */
	case WSATRY_AGAIN:	return "Non-Authoritative: Host not found or SERVERFAIL";	/* 11002 */
	case WSANO_RECOVERY:	return "Non-Recoverable errors, FORMERR, REFUSED, NOTIMP";	/* 11003 */
	case WSANO_DATA:	return "Valid name, no data record of requested type";		/* 11004 */
	case WSAENOMORE:		return "10102: No more results";			/* 10102 */
	case WSAECANCELLED:		return "10103: Call has been canceled";			/* 10103 */
	case WSAEINVALIDPROCTABLE:	return "Procedure call table is invalid";		/* 10104 */
	case WSAEINVALIDPROVIDER:	return "Service provider is invalid";			/* 10105 */
	case WSAEPROVIDERFAILEDINIT:	return "Service provider failed to initialize";		/* 10106 */
	case WSASYSCALLFAILURE:		return "System call failure";				/* 10107 */
	case WSASERVICE_NOT_FOUND:	return "Service not found";				/* 10108 */
	case WSATYPE_NOT_FOUND:		return "Class type not found";				/* 10109 */
	case WSA_E_NO_MORE:		return "10110: No more results";			/* 10110 */
	case WSA_E_CANCELLED:		return "10111: Call was canceled";			/* 10111 */
	case WSAEREFUSED:		return "Database query was refused";			/* 10112 */
	default:
		static s8 _err_unknown[64];
		sprintf(_err_unknown, "Unknown WSAE error (%d)", err);
		return  _err_unknown;
	}
}

static void WINS_GetLocalAddress()
{
	struct hostent *local = NULL;
	s8 buff[MAXHOSTNAMELEN];
	in_addr_t addr;
	s32 err;
	if (myAddr != INADDR_ANY)
		return;
	if (gethostname(buff, MAXHOSTNAMELEN) == SOCKET_ERROR) {
		err = SOCKETERRNO;
		Con_SafePrintf
		    ("WINS_GetLocalAddress: WARNING: gethostname failed (%s)\n",
		     socketerror(err));
		return;
	}
	buff[MAXHOSTNAMELEN - 1] = 0;
	local = gethostbyname(buff);
	err = WSAGetLastError();
	if (local == NULL) {
		Con_SafePrintf
		    ("WINS_GetLocalAddress: gethostbyname failed (%s)\n",
		     __WSAE_StrError(err));
		return;
	}
	myAddr = *(in_addr_t *) local->h_addr_list[0];
	addr = ntohl(myAddr);
	sprintf(my_tcpip_address, "%ld.%ld.%ld.%ld", (addr >> 24) & 0xff,
		(addr >> 16) & 0xff, (addr >> 8) & 0xff, addr & 0xff);
}

sys_socket_t WINS_Init()
{
	s32 i, err;
	s8 buff[MAXHOSTNAMELEN];
	if (COM_CheckParm("-noudp"))
		return -1;
	if (winsock_initialized == 0) {
		err = WSAStartup(MAKEWORD(1, 1), &winsockdata);
		if (err != 0) {
			Con_SafePrintf("Winsock initialization failed (%s)\n",
				       socketerror(err));
			return INVALID_SOCKET;
		}
	}
	winsock_initialized++;
	// determine my name & address
	if (gethostname(buff, MAXHOSTNAMELEN) != 0) {
		err = SOCKETERRNO;
		Con_SafePrintf("WINS_Init: WARNING: gethostname failed (%s)\n",
			       socketerror(err));
	} else {
		buff[MAXHOSTNAMELEN - 1] = 0;
	}
	i = COM_CheckParm("-ip");
	if (i) {
		if (i < com_argc - 1) {
			myAddr = inet_addr(com_argv[i + 1]);
			if (myAddr == INADDR_NONE)
				Sys_Error("%s is not a valid IP address",
					  com_argv[i + 1]);
			strcpy(my_tcpip_address, com_argv[i + 1]);
		} else {
			Sys_Error
			    ("NET_Init: you must specify an IP address after -ip");
		}
	} else {
		myAddr = INADDR_ANY;
		strcpy(my_tcpip_address, "INADDR_ANY");
	}
	if ((net_controlsocket = WINS_OpenSocket(0)) == INVALID_SOCKET) {
		Con_SafePrintf
		    ("WINS_Init: Unable to open control socket, UDP disabled\n");
		if (--winsock_initialized == 0)
			WSACleanup();
		return INVALID_SOCKET;
	}
	broadcastaddr.sin_family = AF_INET;
	broadcastaddr.sin_addr.s_addr = INADDR_BROADCAST;
	broadcastaddr.sin_port = htons((u16)net_hostport);
	Con_SafePrintf("UDP Initialized\n");
	tcpipAvailable = 1;
	return net_controlsocket;
}

void WINS_Shutdown()
{
	WINS_Listen(0);
	WINS_CloseSocket(net_controlsocket);
	if (--winsock_initialized == 0)
		WSACleanup();
}

void WINS_Listen(bool state)
{
	// enable listening
	if (state) {
		if (net_acceptsocket != INVALID_SOCKET)
			return;
		WINS_GetLocalAddress();
		if ((net_acceptsocket =
		     WINS_OpenSocket(net_hostport)) == INVALID_SOCKET)
			Sys_Error("WINS_Listen: Unable to open accept socket");
		return;
	}
	// disable listening
	if (net_acceptsocket == INVALID_SOCKET)
		return;
	WINS_CloseSocket(net_acceptsocket);
	net_acceptsocket = INVALID_SOCKET;
}

sys_socket_t WINS_OpenSocket(s32 port)
{
	sys_socket_t newsocket;
	struct sockaddr_in address;
	u_long _true = 1;
	s32 err;
	if ((newsocket =
	     socket(PF_INET, SOCK_DGRAM, IPPROTO_UDP)) == INVALID_SOCKET) {
		err = SOCKETERRNO;
		Con_SafePrintf("WINS_OpenSocket: %s\n", socketerror(err));
		return INVALID_SOCKET;
	}
	if (ioctlsocket(newsocket, FIONBIO, &_true) == SOCKET_ERROR)
		goto ErrorReturn;
	memset(&address, 0, sizeof(struct sockaddr_in));
	address.sin_family = AF_INET;
	address.sin_addr.s_addr = myAddr;
	address.sin_port = htons((u16)port);
	if (bind(newsocket, (struct sockaddr *)&address, sizeof(address)) == 0)
		return newsocket;
	if (tcpipAvailable) {
		err = SOCKETERRNO;
		Sys_Error("Unable to bind to %s (%s)",
			  WINS_AddrToString((struct qsockaddr *)&address),
			  socketerror(err));
		return INVALID_SOCKET;	/* not reached */
	}
	/* else: we are still in init phase, no need to error */
ErrorReturn:
	err = SOCKETERRNO;
	Con_SafePrintf("WINS_OpenSocket: %s\n", socketerror(err));
	closesocket(newsocket);
	return INVALID_SOCKET;
}

s32 WINS_CloseSocket(sys_socket_t socketid)
{
	if (socketid == net_broadcastsocket)
		net_broadcastsocket = 0;
	return closesocket(socketid);
}
// PartialIPAddress
// this lets you type only as much of the net address as required, using
// the local network components to fill in the rest
static s32 PartialIPAddress(const s8 *in, struct qsockaddr *hostaddr)
{
	s8 buff[256];
	s8 *b;
	s32 addr, mask, num, port, run;
	buff[0] = '.';
	b = buff;
	strcpy(buff + 1, in);
	if (buff[1] == '.')
		b++;
	addr = 0;
	mask = -1;
	while (*b == '.') {
		b++;
		num = 0;
		run = 0;
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
	if (*b++ == ':')
		port = Q_atoi(b);
	else
		port = net_hostport;
	hostaddr->qsa_family = AF_INET;
	((struct sockaddr_in *)hostaddr)->sin_port =
	    htons((u16)port);
	((struct sockaddr_in *)hostaddr)->sin_addr.s_addr =
	    (myAddr & htonl(mask)) | htonl(addr);
	return 0;
}

s32 WINS_Connect(sys_socket_t socketid, struct qsockaddr *addr)
{
	return 0;
}

sys_socket_t WINS_CheckNewConnections()
{
	s8 buf[4096];
	if (net_acceptsocket == INVALID_SOCKET)
		return INVALID_SOCKET;
	if (recvfrom(net_acceptsocket, buf, sizeof(buf), MSG_PEEK, NULL, NULL)
	    != SOCKET_ERROR) {
		return net_acceptsocket;
	}
	return INVALID_SOCKET;
}

s32 WINS_Read(sys_socket_t socketid, u8 *buf, s32 len, struct qsockaddr *addr)
{
	socklen_t addrlen = sizeof(struct qsockaddr);
	s32 ret;
	ret =
	    recvfrom(socketid, (s8 *)buf, len, 0, (struct sockaddr *)addr,
		     &addrlen);
	if (ret == SOCKET_ERROR) {
		s32 err = SOCKETERRNO;
		if (err == NET_EWOULDBLOCK || err == NET_ECONNREFUSED)
			return 0;
		Con_SafePrintf("WINS_Read, recvfrom: %s\n", socketerror(err));
	}
	return ret;
}

static s32 WINS_MakeSocketBroadcastCapable(sys_socket_t socketid)
{
	s32 i = 1;
	// make this socket broadcast capable
	if (setsockopt
	    (socketid, SOL_SOCKET, SO_BROADCAST, (s8 *)&i, sizeof(i))
	    == SOCKET_ERROR) {
		s32 err = SOCKETERRNO;
		Con_SafePrintf("UDP, setsockopt: %s\n", socketerror(err));
		return -1;
	}
	net_broadcastsocket = socketid;
	return 0;
}

s32 WINS_Broadcast(sys_socket_t socketid, u8 *buf, s32 len)
{
	s32 ret;
	if (socketid != net_broadcastsocket) {
		if (net_broadcastsocket != 0)
			Sys_Error
			    ("Attempted to use multiple broadcasts sockets");
		WINS_GetLocalAddress();
		ret = WINS_MakeSocketBroadcastCapable(socketid);
		if (ret == -1) {
			Con_Printf("Unable to make socket broadcast capable\n");
			return ret;
		}
	}
	return WINS_Write(socketid, buf, len,
			  (struct qsockaddr *)&broadcastaddr);
}

s32 WINS_Write(sys_socket_t socketid, u8 *buf, s32 len,struct qsockaddr *addr)
{
	s32 ret;
	ret = sendto(socketid, (s8 *)buf, len, 0, (struct sockaddr *)addr,
		     sizeof(struct qsockaddr));
	if (ret == SOCKET_ERROR) {
		s32 err = SOCKETERRNO;
		if (err == NET_EWOULDBLOCK)
			return 0;
		Con_SafePrintf("WINS_Write, sendto: %s\n", socketerror(err));
	}
	return ret;
}

const s8 *WINS_AddrToString(struct qsockaddr *addr)
{
	static s8 buffer[22];
	s32 haddr;
	haddr = ntohl(((struct sockaddr_in *)addr)->sin_addr.s_addr);
	sprintf(buffer, "%d.%d.%d.%d:%d", (haddr >> 24) & 0xff,
		(haddr >> 16) & 0xff, (haddr >> 8) & 0xff, haddr & 0xff,
		ntohs(((struct sockaddr_in *)addr)->sin_port));
	return buffer;
}

s32 WINS_StringToAddr(const s8 *string, struct qsockaddr *addr)
{
	s32 ha1, ha2, ha3, ha4, hp, ipaddr;
	sscanf(string, "%d.%d.%d.%d:%d", &ha1, &ha2, &ha3, &ha4, &hp);
	ipaddr = (ha1 << 24) | (ha2 << 16) | (ha3 << 8) | ha4;
	addr->qsa_family = AF_INET;
	((struct sockaddr_in *)addr)->sin_addr.s_addr = htonl(ipaddr);
	((struct sockaddr_in *)addr)->sin_port = htons((u16)hp);
	return 0;
}

s32 WINS_GetSocketAddr(sys_socket_t socketid, struct qsockaddr *addr)
{
	socklen_t addrlen = sizeof(struct qsockaddr);
	in_addr_t a;
	memset(addr, 0, sizeof(struct qsockaddr));
	getsockname(socketid, (struct sockaddr *)addr, &addrlen);
	a = ((struct sockaddr_in *)addr)->sin_addr.s_addr;
	if (a == 0 || a == htonl(INADDR_LOOPBACK))
		((struct sockaddr_in *)addr)->sin_addr.s_addr = myAddr;
	return 0;
}

s32 WINS_GetNameFromAddr(struct qsockaddr *addr, s8 *name)
{
	struct hostent *hostentry;
	hostentry =
	    gethostbyaddr((s8 *)&((struct sockaddr_in *)addr)->sin_addr,
			  sizeof(struct in_addr), AF_INET);
	if (hostentry) {
		Q_strncpy(name, (s8 *)hostentry->h_name, NET_NAMELEN - 1);
		return 0;
	}
	Q_strcpy(name, WINS_AddrToString(addr));
	return 0;
}

s32 WINS_GetAddrFromName(const s8 *name, struct qsockaddr *addr)
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

s32 WINS_AddrCompare(struct qsockaddr *addr1, struct qsockaddr *addr2)
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

s32 WINS_GetSocketPort(struct qsockaddr *addr)
{
	return ntohs(((struct sockaddr_in *)addr)->sin_port);
}

s32 WINS_SetSocketPort(struct qsockaddr *addr, s32 port)
{
	((struct sockaddr_in *)addr)->sin_port = htons((u16)port);
	return 0;
}
