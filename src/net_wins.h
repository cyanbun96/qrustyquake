// Copyright (C) 1996-2001 Id Software, Inc.
// Copyright (C) 2010-2014 QuakeSpasm developers
// GPLv3 See LICENSE for details.

sys_socket_t WINS_Init();
void WINS_Shutdown();
void WINS_Listen(bool state);
sys_socket_t WINS_OpenSocket(s32 port);
s32 WINS_CloseSocket(sys_socket_t socketid);
s32 WINS_Connect(sys_socket_t socketid, struct qsockaddr *addr);
sys_socket_t WINS_CheckNewConnections();
s32 WINS_Read(sys_socket_t socketid, u8 *buf, s32 len, struct qsockaddr *addr);
s32 WINS_Write(sys_socket_t socketid, u8 *buf, s32 len, struct qsockaddr *addr);
s32 WINS_Broadcast(sys_socket_t socketid, u8 *buf, s32 len);
const s8 *WINS_AddrToString(struct qsockaddr *addr);
s32 WINS_StringToAddr(const s8 *string, struct qsockaddr *addr);
s32 WINS_GetSocketAddr(sys_socket_t socketid, struct qsockaddr *addr);
s32 WINS_GetNameFromAddr(struct qsockaddr *addr, s8 *name);
s32 WINS_GetAddrFromName(const s8 *name, struct qsockaddr *addr);
s32 WINS_AddrCompare(struct qsockaddr *addr1, struct qsockaddr *addr2);
s32 WINS_GetSocketPort(struct qsockaddr *addr);
s32 WINS_SetSocketPort(struct qsockaddr *addr, s32 port);
