package stego;

import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.charset.CodingErrorAction;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.zip.CRC32;

/** Streaming file manifest and disk-backed, order-independent assembly.
 * @author Yaoyao */
public final class StegoArchive{

	private StegoArchive(){}
	public static final int MAX_NAME_BYTES=65535;

	public static void validateName(String name){
		if(name==null || name.isEmpty() || name.equals(".") || name.equals("..") ||
			name.indexOf('/')>=0 || name.indexOf('\\')>=0 || name.indexOf(':')>=0 || name.indexOf('\0')>=0){
			throw new IllegalArgumentException("Metadata filename must be a plain basename: "+name);
		}
	}

	/** Count parts by base-4 part-number width, without walking every future read. */
	public static long countParts(long bytes, int length, long group){
		long part=0;
		while(bytes>0){
			int cap=StegoCodec.capacity(length, group, part);
			if(cap<1){throw new IllegalArgumentException("Read length too short for metadata envelope");}
			int digits=StegoCodec.integerBases(part)-3;
			long boundary=digits>=32 ? Long.MAX_VALUE : 1L<<(2*digits);
			long span=boundary-part, needed=(bytes-1)/cap+1;
			if(span==0){throw new IllegalArgumentException("Attachment requires too many parts");}
			if(needed<=span){return Math.addExact(part, needed);}
			bytes-=span*cap;
			part=boundary;
		}
		return part;
	}

	/** File stream: u32 name bytes, UTF-8 basename, u64 parts, u64 bytes, data, CRC32. */
	public static final class Payload extends InputStream{
		public Payload(String name, long size, InputStream input, int length, long group) throws IOException{
			validateName(name);
			byte[] nameBytes=name.getBytes(StandardCharsets.UTF_8);
			if(nameBytes.length>MAX_NAME_BYTES || size<0){throw new IllegalArgumentException("Invalid filename or file size");}
			parts=countParts(Math.addExact(size, 24L+nameBytes.length), length, group);
			ByteArrayOutputStream buffer=new ByteArrayOutputStream();
			DataOutputStream out=new DataOutputStream(buffer);
			out.writeInt(nameBytes.length);
			out.write(nameBytes);
			out.writeLong(parts);
			out.writeLong(size);
			header=buffer.toByteArray();
			raw=input;
			remaining=size;
		}
		@Override
		public int read() throws IOException{
			int n=read(single, 0, 1);
			return n<0 ? -1 : single[0]&255;
		}
		@Override
		public int read(byte[] b, int off, int len) throws IOException{
			if(len==0){return 0;}
			if(headerPos<header.length){
				int n=Math.min(len, header.length-headerPos);
				System.arraycopy(header, headerPos, b, off, n);
				headerPos+=n;
				return n;
			}
			if(remaining>0){
				int n=raw.read(b, off, (int)Math.min(len, remaining));
				if(n<0){throw new EOFException("Input file shrank while encoding");}
				crc.update(b, off, n);
				remaining-=n;
				return n;
			}
			if(trailer==null){
				if(raw.read()!=-1){throw new IOException("Input file grew while encoding");}
				long value=crc.getValue();
				trailer=new byte[]{(byte)(value>>>24),(byte)(value>>>16),(byte)(value>>>8),(byte)value};
			}
			if(trailerPos==4){return -1;}
			int n=Math.min(len, 4-trailerPos);
			System.arraycopy(trailer, trailerPos, b, off, n);
			trailerPos+=n;
			return n;
		}
		@Override
		public void close() throws IOException{raw.close();}
		public final long parts;
		private final InputStream raw;
		private final byte[] header, single=new byte[1];
		private final CRC32 crc=new CRC32();
		private byte[] trailer;
		private int headerPos, trailerPos;
		private long remaining;
	}

	public static final class Decoder implements AutoCloseable{
		public Decoder(Path output, int maxFiles) throws IOException{
			Files.createDirectories(output);
			directory=output.toRealPath();
			work=Files.createTempDirectory(directory, ".stegofast-");
			limit=maxFiles;
		}
		public void accept(StegoCodec.Chunk chunk) throws IOException{
			Group group=last;
			if(group==null || group.id!=chunk.group){
				group=null;
				for(Group g : groups){if(g.id==chunk.group){group=g; break;}}
				if(group==null){
					if(groups.size()>=limit){throw new IOException("Too many metadata files; increase maxfiles");}
					group=new Group(chunk.group, chunk.length, work);
					groups.add(group);
				}
				last=group;
			}
			group.accept(chunk);
		}
		/** Validate every group and output name before publishing any file. Never overwrites. */
		public int finish() throws IOException{
			HashSet<String> names=new HashSet<String>();
			for(Group g : groups){
				g.prepare();
				if(!names.add(g.name) || Files.exists(directory.resolve(g.name), java.nio.file.LinkOption.NOFOLLOW_LINKS)){
					throw new IOException("Duplicate or existing extraction filename: "+g.name);
				}
			}
			for(Group g : groups){Files.move(g.product, directory.resolve(g.name));}
			close();
			// Spools are regenerable; the validated recovered files now retain their data.
			for(Group g : groups){Files.delete(g.spool);}
			Files.delete(work);
			return groups.size();
		}
		@Override
		public void close() throws IOException{
			IOException failure=null;
			for(Group g : groups){try{g.file.close();}catch(IOException e){failure=e;}}
			if(failure!=null){throw failure;}
		}
		public Path workDirectory(){return work;}
		private final Path directory, work;
		private final int limit;
		private final ArrayList<Group> groups=new ArrayList<Group>();
		private Group last;
	}

	private static final class Group{
		Group(long group, int readLength, Path work) throws IOException{
			id=group; length=readLength; stride=4+readLength/4;
			spool=work.resolve(Long.toHexString(id)+".parts");
			product=work.resolve(Long.toHexString(id)+".decoded");
			file=new RandomAccessFile(spool.toFile(), "rw");
		}
		void accept(StegoCodec.Chunk chunk) throws IOException{
			if(chunk.length!=length){throw new IOException("Conflicting read lengths for metadata group "+id);}
			long offset=Math.multiplyExact(chunk.part, (long)stride);
			file.seek(offset);
			int stored=Math.addExact(offset, 4)<=file.length() ? file.readInt() : 0;
			if(stored!=0){
				if(stored!=chunk.data.length){throw new IOException("Conflicting duplicate metadata part");}
				byte[] previous=new byte[stored];
				file.readFully(previous);
				if(!Arrays.equals(previous, chunk.data)){throw new IOException("Conflicting duplicate metadata part");}
				return;
			}
			file.seek(offset);
			file.writeInt(chunk.data.length);
			file.write(chunk.data);
			count++;
			maxPart=Math.max(maxPart, chunk.part);
		}
		void prepare() throws IOException{
			if(count!=Math.addExact(maxPart, 1)){throw new IOException("Missing parts in metadata group "+id);}
			try(DataInputStream in=new DataInputStream(new PartInput(this))){
				int size=in.readInt();
				if(size<=0 || size>MAX_NAME_BYTES){throw new IOException("Invalid metadata filename length");}
				byte[] bytes=new byte[size];
				in.readFully(bytes);
				name=StandardCharsets.UTF_8.newDecoder().onMalformedInput(CodingErrorAction.REPORT)
					.onUnmappableCharacter(CodingErrorAction.REPORT).decode(ByteBuffer.wrap(bytes)).toString();
				validateName(name);
				long parts=in.readLong(), remaining=in.readLong();
				if(parts!=count || remaining<0){throw new IOException("Metadata manifest count mismatch");}
				CRC32 crc=new CRC32();
				try(OutputStream out=Files.newOutputStream(product, StandardOpenOption.CREATE_NEW)){
					byte[] buffer=new byte[65536];
					while(remaining>0){
						int n=in.read(buffer, 0, (int)Math.min(buffer.length, remaining));
						if(n<0){throw new EOFException("Truncated metadata file");}
						out.write(buffer, 0, n);
						crc.update(buffer, 0, n);
						remaining-=n;
					}
				}
				if((in.readInt()&0xffffffffL)!=crc.getValue() || in.read()!=-1){
					throw new IOException("Metadata file CRC or byte count mismatch");
				}
			}
		}
		final long id;
		final int length, stride;
		final Path spool, product;
		final RandomAccessFile file;
		long count, maxPart=-1;
		String name;
	}

	private static final class PartInput extends InputStream{
		PartInput(Group g){group=g; buffer=new byte[g.stride-4];}
		@Override
		public int read() throws IOException{
			if(!availablePart()){return -1;}
			return buffer[pos++]&255;
		}
		@Override
		public int read(byte[] b, int off, int len) throws IOException{
			if(len==0){return 0;}
			if(!availablePart()){return -1;}
			int n=Math.min(len, size-pos);
			System.arraycopy(buffer, pos, b, off, n);
			pos+=n;
			return n;
		}
		private boolean availablePart() throws IOException{
			if(pos<size){return true;}
			if(part>=group.count){return false;}
			group.file.seek(Math.multiplyExact(part++, (long)group.stride));
			size=group.file.readInt();
			if(size<=0 || size>buffer.length){throw new IOException("Missing or invalid spooled part");}
			group.file.readFully(buffer, 0, size);
			pos=0;
			return true;
		}
		final Group group;
		final byte[] buffer;
		int pos, size;
		long part;
	}
}
