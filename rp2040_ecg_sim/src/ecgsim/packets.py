"""Common binary packet format used by ALL architectures (same sample representation).

Layout (little-endian), 10 + 2*N bytes:
    0  u16  magic 0x5AA5
    2  u32  sequence number of first sample slot
    6  u8   number of sample slots N
    7  u8   flags (bit0: packet contains acquisition-missed slots)
    8  u16[N] samples: 12-bit ADC code, or 0xFFFF = "no sample acquired for this slot"
    8+2N u16 CRC-16/CCITT-FALSE (poly 0x1021, init 0xFFFF) over bytes [0, 8+2N)
CRC-16-CCITT is chosen because the RP2040 DMA sniffer can compute it in hardware.
"""
from __future__ import annotations

import binascii
import struct

import numpy as np

MAGIC = 0x5AA5
MISSING = 0xFFFF
HDR = struct.Struct("<HIBB")


def crc16(data: bytes) -> int:
    return binascii.crc_hqx(data, 0xFFFF)


def build_packet(seq0: int, codes: np.ndarray) -> bytes:
    codes = np.asarray(codes, dtype="<u2")
    flags = 1 if np.any(codes == MISSING) else 0
    body = HDR.pack(MAGIC, seq0 & 0xFFFFFFFF, len(codes), flags) + codes.tobytes()
    return body + struct.pack("<H", crc16(body))


def parse_packet(pkt: bytes):
    """Return (seq0, codes) or None if the packet fails framing/CRC checks."""
    if len(pkt) < HDR.size + 2:
        return None
    body, (crc,) = pkt[:-2], struct.unpack("<H", pkt[-2:])
    if crc16(body) != crc:
        return None
    magic, seq0, n, _flags = HDR.unpack(body[: HDR.size])
    if magic != MAGIC or len(body) != HDR.size + 2 * n:
        return None
    return seq0, np.frombuffer(body[HDR.size:], dtype="<u2")


def packet_len(n: int) -> int:
    return HDR.size + 2 * n + 2


def flip_bits(pkt: bytes, positions) -> bytes:
    b = bytearray(pkt)
    for p in positions:
        b[p // 8] ^= 1 << (p % 8)
    return bytes(b)
