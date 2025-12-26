// Copyright (C) 1996-2001 Id Software, Inc.
// Copyright (C) 2010-2014 QuakeSpasm developers
// GPLv3 See LICENSE for details.

#include "quakedef.h"

static bool localconnectpending = 0;
static qsocket_t *loop_client = NULL;
static qsocket_t *loop_server = NULL;

s32 Loop_Init()
{
	if (cls.state == ca_dedicated)
		return -1;
	return 0;
}

void Loop_Shutdown()
{
}

void Loop_Listen(bool state)
{
	(void)state; // removing the argument would break things, probably
} // mentioning it here with void just to quiet the compiler

void Loop_SearchForHosts(bool xmit)
{
	(void)xmit; // same as above
	if (!sv.active)
		return;
	hostCacheCount = 1;
	if (Q_strcmp(hostname.string, "UNNAMED") == 0)
		Q_strcpy(hostcache[0].name, "local");
	else
		Q_strcpy(hostcache[0].name, hostname.string);
	Q_strcpy(hostcache[0].map, sv.name);
	hostcache[0].users = net_activeconnections;
	hostcache[0].maxusers = svs.maxclients;
	hostcache[0].driver = net_driverlevel;
	Q_strcpy(hostcache[0].cname, "local");
}

qsocket_t *Loop_Connect(const s8 *host)
{
	if (Q_strcmp(host, "local") != 0)
		return NULL;
	localconnectpending = 1;
	if (!loop_client) {
		if ((loop_client = NET_NewQSocket()) == NULL) {
			Con_Printf("Loop_Connect: no qsocket available\n");
			return NULL;
		}
		Q_strcpy(loop_client->address, "localhost");
	}
	loop_client->receiveMessageLength = 0;
	loop_client->sendMessageLength = 0;
	loop_client->canSend = 1;
	if (!loop_server) {
		if ((loop_server = NET_NewQSocket()) == NULL) {
			Con_Printf("Loop_Connect: no qsocket available\n");
			return NULL;
		}
		Q_strcpy(loop_server->address, "LOCAL");
	}
	loop_server->receiveMessageLength = 0;
	loop_server->sendMessageLength = 0;
	loop_server->canSend = 1;
	loop_client->driverdata = (void *)loop_server;
	loop_server->driverdata = (void *)loop_client;
	return loop_client;
}

qsocket_t *Loop_CheckNewConnections()
{
	if (!localconnectpending)
		return NULL;
	localconnectpending = 0;
	loop_server->sendMessageLength = 0;
	loop_server->receiveMessageLength = 0;
	loop_server->canSend = 1;
	loop_client->sendMessageLength = 0;
	loop_client->receiveMessageLength = 0;
	loop_client->canSend = 1;
	return loop_server;
}

static s32 IntAlign(s32 value)
{
	return (value + (sizeof(s32) - 1)) & (~(sizeof(s32) - 1));
}

s32 Loop_GetMessage(qsocket_t *sock)
{
	if (sock->receiveMessageLength == 0)
		return 0;
	s32 ret = sock->receiveMessage[0];
	s32 length = sock->receiveMessage[1] + (sock->receiveMessage[2] << 8);
	// alignment u8 skipped here
	SZ_Clear(&net_message);
	SZ_Write(&net_message, &sock->receiveMessage[4], length);
	length = IntAlign(length + 4);
	sock->receiveMessageLength -= length;
	if (sock->receiveMessageLength)
		memmove(sock->receiveMessage, &sock->receiveMessage[length],
			sock->receiveMessageLength);
	if (sock->driverdata && ret == 1)
		((qsocket_t *) sock->driverdata)->canSend = 1;
	return ret;
}

s32 Loop_SendMessage(qsocket_t *sock, sizebuf_t *data)
{
	if (!sock->driverdata)
		return -1;
	s32 *bufferLength = &((qsocket_t *) sock->driverdata)->receiveMessageLength;
	if ((*bufferLength + data->cursize + 4) > NET_MAXMESSAGE)
		Sys_Error("Loop_SendMessage: overflow");
	u8 *buffer =
	    ((qsocket_t *) sock->driverdata)->receiveMessage + *bufferLength;
	*buffer++ = 1; // message type
	*buffer++ = data->cursize & 0xff; // length
	*buffer++ = data->cursize >> 8;
	buffer++; // align
	Q_memcpy(buffer, data->data, data->cursize); // message
	*bufferLength = IntAlign(*bufferLength + data->cursize + 4);
	sock->canSend = 0;
	return 1;
}

s32 Loop_SendUnreliableMessage(qsocket_t *sock, sizebuf_t *data)
{
	if (!sock->driverdata)
		return -1;
	s32 *bufferLength =
		&((qsocket_t *) sock->driverdata)->receiveMessageLength;
	if ((*bufferLength + data->cursize + sizeof(u8) + sizeof(s16)) >
	    NET_MAXMESSAGE)
		return 0;
	u8 *buffer =
	    ((qsocket_t *) sock->driverdata)->receiveMessage + *bufferLength;
	*buffer++ = 2; // message type
	*buffer++ = data->cursize & 0xff; // length
	*buffer++ = data->cursize >> 8;
	buffer++; // align
	Q_memcpy(buffer, data->data, data->cursize); // message
	*bufferLength = IntAlign(*bufferLength + data->cursize + 4);
	return 1;
}

bool Loop_CanSendMessage(qsocket_t *sock)
{
	if (!sock->driverdata)
		return 0;
	return sock->canSend;
}

bool Loop_CanSendUnreliableMessage(qsocket_t *sock)
{
	(void)sock;
	return 1;
}

void Loop_Close(qsocket_t *sock)
{
	if (sock->driverdata)
		((qsocket_t *) sock->driverdata)->driverdata = NULL;
	sock->receiveMessageLength = 0;
	sock->sendMessageLength = 0;
	sock->canSend = 1;
	if (sock == loop_client)
		loop_client = NULL;
	else
		loop_server = NULL;
}
