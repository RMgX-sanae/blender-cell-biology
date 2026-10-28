"""Validate cached Blender PNGs before a render resumes; standard library only."""
import struct
import zlib


def valid_png(path, width, height):
    try:
        data=path.read_bytes()
        if data[:8]!=b'\x89PNG\r\n\x1a\n':return False
        position=8;compressed=[];expected=None;ended=False
        while position+12<=len(data):
            size=struct.unpack('>I',data[position:position+4])[0]
            kind=data[position+4:position+8];end=position+8+size
            if end+4>len(data):return False
            chunk=data[position+8:end]
            if zlib.crc32(kind+chunk)&0xffffffff!=struct.unpack('>I',data[end:end+4])[0]:return False
            if kind==b'IHDR':
                w,h,bits,color,compression,filtering,interlace=struct.unpack('>IIBBBBB',chunk)
                channels={0:1,2:3,4:2,6:4}.get(color)
                if (w,h)!=(width,height) or not channels or bits not in {8,16} or interlace:return False
                expected=height*(1+width*channels*bits//8)
            elif kind==b'IDAT':compressed.append(chunk)
            elif kind==b'IEND':ended=True;break
            position=end+4
        if not ended or expected is None:return False
        decoder=zlib.decompressobj();raw=decoder.decompress(b''.join(compressed),expected+1)
        return len(raw)==expected and decoder.eof
    except (OSError,ValueError,struct.error,zlib.error):
        return False
