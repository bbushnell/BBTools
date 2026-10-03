package prot;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.ByteBuffer;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.security.DigestInputStream;
import java.security.DigestOutputStream;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.Arrays;

import fileIO.ReadWrite;
import parse.Parser;
import stream.bam.BgzfOutputStream;
import structures.ByteBuilder;

/** Lossless, self-described A48 transport for a native MQHB byte image.
 * Only nonnegative payload integers change representation. Metadata, float
 * settings, consensus digests and native integrity trailer remain byte-exact.
 * @author Yoimiya
 */
public final class HbmA48Codec {

	/** Converts a fresh output with mode=encode|decode; no existing model is overwritten. */
	public static void main(String[] args) throws IOException{
		if(args.length==1 && args[0].equalsIgnoreCase("selftest=t")){HbmA48CodecTest.main(new String[0]); return;}
		String in=null, out=null, mode=null;
		for(String arg:Parser.parseConfig(args)){
			final int split=arg.indexOf('=');
			if(split<1 || split==arg.length()-1){throw new IllegalArgumentException("Expected nonempty flag=value");}
			final String key=arg.substring(0, split), value=arg.substring(split+1);
			if(key.equalsIgnoreCase("in") && in==null){in=value;}
			else if(key.equalsIgnoreCase("out") && out==null){out=value;}
			else if(key.equalsIgnoreCase("mode") && mode==null){mode=value;}
			else{throw new IllegalArgumentException("Unknown or duplicate option: "+key);}
		}
		if(in==null || out==null || !("encode".equals(mode) || "decode".equals(mode))){
			throw new IllegalArgumentException("Require in= out=fresh-file mode=encode|decode");
		}
		if(mode.equals("encode")){encode(Paths.get(in), Paths.get(out));}
		else{decode(Paths.get(in), Paths.get(out));}
		System.out.println("HBM_A48_PASS mode="+mode+" bytes="+Files.size(Paths.get(out)));
	}

	private HbmA48Codec(){}

	/** Detects the gzip transport only; the decoder still validates the critical header. */
	static boolean isEncoded(final Path path) throws IOException{
		try(InputStream in=Files.newInputStream(path)){return in.read()==31 && in.read()==139;}
	}

	/** Writes BGZF level6 with an explicit A48 header and one canonical EOF marker. */
	public static void encode(final Path input, final Path output) throws IOException{
		fresh(output);
		final MessageDigest digest=sha256();
		try(DataInputStream in=new DataInputStream(new DigestInputStream(
				new BufferedInputStream(Files.newInputStream(input), 1<<20), digest))){
			final OutputStream target=Files.newOutputStream(output, StandardOpenOption.CREATE_NEW);
			boolean complete=false;
			try(OutputStream owned=target;
					BgzfOutputStream compressed=new BgzfOutputStream(new BufferedOutputStream(owned, 1<<20), 6)){
				final DataOutputStream out=new DataOutputStream(compressed);
				out.write(HEADER);
				final long[] lengths=metadata(in, out);
				final ByteBuilder buffer=new ByteBuilder(1<<16);
				final byte[] consensus=new byte[32];
				for(long length:lengths){
					append(buffer, in.readInt());
					in.readFully(consensus); buffer.append(consensus);
					for(long position=36; position<length; position+=4){
						append(buffer, in.readInt());
						if(buffer.length()>=65536){out.write(buffer.array, 0, buffer.length()); buffer.clear();}
					}
					out.write(buffer.array, 0, buffer.length()); buffer.clear();
				}
				final byte[] prefix=digest.digest(), trailer=new byte[36];
				in.readFully(trailer); verifyTrailer(prefix, trailer);
				if(in.read()!=-1){throw new IOException("Trailing MQHB input bytes");}
				out.write(trailer); out.flush(); compressed.writeEOF();
				complete=true;
			}finally{if(!complete){Files.deleteIfExists(output);}}
		}
	}

	/** Restores the native byte image, validating its digest before returning it. */
	public static void decode(final Path input, final Path output) throws IOException{
		fresh(output);
		final MessageDigest digest=sha256();
		// The transport uses .bgz, so request gzip explicitly instead of suffix
		// dispatch. Buffer decoded bytes: integer() consumes one byte at a time.
		// Native BGZF can inflate blocks in parallel; subprocesses are disallowed.
		try(DataInputStream in=new DataInputStream(new BufferedInputStream(
				ReadWrite.getGZipInputStream(input.toString(), false, true), 1<<16))){
			final OutputStream target=Files.newOutputStream(output, StandardOpenOption.CREATE_NEW);
			boolean complete=false;
			try(OutputStream owned=target;
					DataOutputStream out=new DataOutputStream(new BufferedOutputStream(new DigestOutputStream(owned, digest), 1<<20))){
				final byte[] header=new byte[HEADER.length]; in.readFully(header);
				if(!Arrays.equals(header, HEADER)){throw new IOException("Unsupported HBM A48 header/version/encoding");}
				final long[] lengths=metadata(in, out);
				final byte[] consensus=new byte[32];
				for(long length:lengths){
					out.writeInt(integer(in));
					in.readFully(consensus); out.write(consensus);
					for(long position=36; position<length; position+=4){out.writeInt(integer(in));}
				}
				// Flush the complete native prefix through the digest before finalizing
				// it. Buffering above DigestOutputStream batches the four-byte integers.
				out.flush();
				final byte[] prefix=digest.digest(), trailer=new byte[36];
				in.readFully(trailer); verifyTrailer(prefix, trailer);
				if(in.read()!=-1){throw new IOException("Trailing encoded HBM bytes");}
				out.write(trailer); out.flush();
				complete=true;
			}finally{if(!complete){Files.deleteIfExists(output);}}
		}
	}

	/** Copies native metadata, enforcing allocation and reconstructed-size limits first. */
	private static long[] metadata(final DataInputStream in, final DataOutputStream out) throws IOException{
		final byte[] header=new byte[HbmBundleFormat.FIXED_HEADER_LEN]; in.readFully(header);
		final ByteBuffer fields=ByteBuffer.wrap(header);
		if(!HbmBundleFormat.magicMatches(header, 0) || fields.getShort(4)!=1 || fields.getShort(6)!=22 || fields.getShort(8)!=0){
			throw new IOException("A48 transport requires native MQHBv1 with22 residues and pad0");
		}
		final int families=fields.getInt(10);
		if(families<1 || families>100000){throw new IOException("HBM family count exceeds A48 directory limit");}
		out.write(header);
		final long[] lengths=new long[families], offsets=new long[families];
		long position=header.length;
		for(int i=0; i<families; i++){
			final int size=in.readInt();
			if(size<1 || size>HbmBundleFormat.CAP_REPID){throw new IOException("Invalid representative ID length");}
			out.writeInt(size);
			final byte[] entry=new byte[size+20]; in.readFully(entry); out.write(entry);
			final ByteBuffer fieldsEntry=ByteBuffer.wrap(entry); fieldsEntry.position(size);
			offsets[i]=fieldsEntry.getLong(); lengths[i]=fieldsEntry.getLong();
			if(lengths[i]<36 || lengths[i]>HbmBundleLoader.CAP_BLOCK_BYTES || (lengths[i]-36)%4!=0){
				throw new IOException("Invalid native block length at family "+i);
			}
			position+=size+24;
		}
		for(int i=0; i<families; i++){
			if(offsets[i]!=position || lengths[i]>HbmBundleFormat.CAP_TOTAL_BYTES-36-position){
				throw new IOException("Noncontiguous or oversized decoded HBM at family "+i);
			}
			position+=lengths[i];
		}
		return lengths;
	}

	/** Canonical native six-bit alphabet, TAB-delimited; counts are never rounded. */
	private static void append(final ByteBuilder out, final int count) throws IOException{
		if(count<0){throw new IOException("Negative HBM payload integer");}
		out.appendA48(count).tab();
	}

	/** Rejects ambiguous tokens, non-alphabet bytes and signed32-bit overflow. */
	static int integer(final DataInputStream in) throws IOException{
		int value=0, count=0, first=-1;
		for(int symbol=in.readUnsignedByte(); symbol!='\t'; symbol=in.readUnsignedByte()){
			if(symbol<48 || symbol>111 || count==6 || value>(Integer.MAX_VALUE-(symbol-48))/64){
				throw new IOException("Invalid or overflowing A48 HBM integer");
			}
			if(count==0){first=symbol;}
			value=(value<<6)|(symbol-48); count++;
		}
		if(count==0 || (count>1 && first=='0')){throw new IOException("Noncanonical A48 HBM integer");}
		assert(value>=0) : "A48 counts must fit native nonnegative signed32-bit fields";
		return value;
	}

	/** The native trailer protects every reconstructed metadata and payload byte. */
	private static void verifyTrailer(final byte[] prefix, final byte[] trailer) throws IOException{
		for(int i=0; i<prefix.length; i++){
			if(prefix[i]!=trailer[i]){throw new IOException("Native HBM prefix digest mismatch");}
		}
		if(!HbmBundleFormat.magicMatches(trailer, 32)){throw new IOException("Native HBM trailer magic mismatch");}
	}

	/** Fresh output is mandatory for preservation of accepted resources. */
	private static void fresh(final Path output) throws IOException{
		if(Files.exists(output)){throw new IOException("Output already exists: "+output);}
	}
	private static MessageDigest sha256(){
		try{return MessageDigest.getInstance("SHA-256");}
		catch(NoSuchAlgorithmException e){throw new IllegalStateException("Missing SHA-256 implementation", e);}
	}
	private static final byte[] HEADER="MQHA\t1\n#encoding\tA48\n".getBytes(StandardCharsets.US_ASCII);
}
