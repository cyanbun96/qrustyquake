// Copyright (C) 1996-1997 Id Software, Inc. GPLv3 See LICENSE for details.
#include "quakedef.h"
net_driver_t net_drivers[] = {
	{ "Loopback", 0, Loop_Init, Loop_Listen, Loop_SearchForHosts,
	 Loop_Connect, Loop_CheckNewConnections, Loop_GetMessage,
	 Loop_SendMessage, Loop_SendUnreliableMessage, Loop_CanSendMessage,
	 Loop_CanSendUnreliableMessage, Loop_Close, Loop_Shutdown },
	{ "Datagram", 0, Datagram_Init, Datagram_Listen,
	 Datagram_SearchForHosts, Datagram_Connect,
	 Datagram_CheckNewConnections, Datagram_GetMessage,
	 Datagram_SendMessage, Datagram_SendUnreliableMessage,
	 Datagram_CanSendMessage, Datagram_CanSendUnreliableMessage,
	 Datagram_Close, Datagram_Shutdown }
};
const s32 net_numdrivers = Q_COUNTOF(net_drivers);
net_landriver_t net_landrivers[] = {
	{ "UDP", 0, 0, UDP_Init, UDP_Shutdown, UDP_Listen, UDP_OpenSocket,
	 UDP_CloseSocket, UDP_Connect, UDP_CheckNewConnections, UDP_Read,
	 UDP_Write, UDP_Broadcast, UDP_AddrToString, UDP_StringToAddr,
	 UDP_GetSocketAddr, UDP_GetNameFromAddr, UDP_GetAddrFromName,
	 UDP_AddrCompare, UDP_GetSocketPort, UDP_SetSocketPort } };
const s32 net_numlandrivers = Q_COUNTOF(net_landrivers);
