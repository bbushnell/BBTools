package prot;

import java.io.BufferedOutputStream;
import java.io.Closeable;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.security.MessageDigest;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import dna.AminoAcid;
import parse.Parser;
import stream.bam.BgzfOutputStream;
import structures.ByteBuilder;

/**
 * Experimental all-text dump of native MQHB HBM bundles.  The row stream preserves
 * the native REF/DEL insertion topology while writing only retained candidate
 * states and, when supplied, pre-pruning row totals from the corresponding
 * original bundle.
 *
 * @author Yelan
 */
public final class HbmTextDump {

	private HbmTextDump(){}

	private static final String FORMAT="hbm_text_v1";
	private static final int DENSE_SYMBOL_THRESHOLD=(HbmBundleFormat.NAA+1)/2;

	public static void main(final String[] args) throws Exception{
		if(args.length==1 && args[0].equalsIgnoreCase("selftest=t")){
			HbmTextDumpTest.main(new String[0]);
			return;
		}
		String input=null, original=null, reference=null, output=null;
		String encoding="sparse", coords="explicit";
		int minCount=1;
		boolean bgz=false;
		for(final String arg : Parser.parseConfig(args)){
			final int equals=arg.indexOf('=');
			if(equals<1 || equals==arg.length()-1){throw new IllegalArgumentException("Expected flag=value: "+arg);}
			final String key=arg.substring(0, equals), value=arg.substring(equals+1);
			if(key.equalsIgnoreCase("in") && input==null){input=value;}
			else if(key.equalsIgnoreCase("orig") && original==null){original=value;}
			else if(key.equalsIgnoreCase("ref") && reference==null){reference=value;}
			else if(key.equalsIgnoreCase("out") && output==null){output=value;}
			else if(key.equalsIgnoreCase("encoding")){encoding=value;}
			else if(key.equalsIgnoreCase("coords")){coords=value;}
			else if(key.equalsIgnoreCase("mincount")){minCount=parseInt(value, key);}
			else if(key.equalsIgnoreCase("bgz")){bgz=parseBoolean(value, key);}
			else{throw new IllegalArgumentException("Unknown or duplicate option: "+key);}
		}
		if(input==null || reference==null || output==null){
			throw new IllegalArgumentException("Required in= ref= out=fresh-file [orig=] [encoding=sparse|dense|adaptive] [coords=explicit|implicit] [mincount=1] [bgz=t]");
		}
		run(Paths.get(input), original==null ? null : Paths.get(original),
			reference, Paths.get(output), Encoding.parse(encoding), CoordMode.parse(coords),
			minCount, bgz);
	}

	static void run(final Path input, final Path original, final String reference,
			final Path output, final Encoding encoding, final CoordMode coordMode,
			final int minCount, final boolean bgz) throws Exception{
		if(minCount<1){throw new IllegalArgumentException("mincount must be at least 1: "+minCount);}
		final HashMap<String, byte[]> consensus=loadConsensus(reference);
		final Bundle candidate=Bundle.load(input, consensus);
		final Bundle unpruned=(original==null ? candidate : Bundle.load(original, consensus));
		candidate.assertCompatible(unpruned);
		write(candidate, unpruned, output, encoding, coordMode, minCount, bgz);
	}

	private static HashMap<String, byte[]> loadConsensus(final String reference) throws IOException{
		final List<ProteinSequence> sequences=ProteinSearch.readFasta(reference);
		if(sequences.isEmpty()){throw new IllegalArgumentException("Empty consensus reference: "+reference);}
		final HashMap<String, byte[]> map=new HashMap<String, byte[]>(sequences.size()*2);
		for(final ProteinSequence sequence : sequences){
			if(map.put(sequence.id, sequence.enc)!=null){throw new IllegalArgumentException("Duplicate consensus: "+sequence.id);}
		}
		return map;
	}

	private static void write(final Bundle candidate, final Bundle original, final Path output,
			final Encoding encoding, final CoordMode coordMode,
			final int minCount, final boolean bgz) throws Exception{
		fresh(output);
		final OutputStream raw=Files.newOutputStream(output, StandardOpenOption.CREATE_NEW);
		boolean complete=false;
		try{
			final OutputStream owned=raw;
			final OutputStream compressed=(bgz ?
				new BgzfOutputStream(new BufferedOutputStream(owned, 1<<20), 6) :
				new BufferedOutputStream(owned, 1<<20));
			try(TextOut out=new TextOut(compressed)){
				out.append("#format\t").append(FORMAT).nl();
				out.append("#candidate_sha80\t").append(candidate.sourceSha80).nl();
				out.append("#original_sha80\t").append(original.sourceSha80).nl();
				out.append("#sum_source\t").append(candidate==original ? "candidate" : "original").nl();
				out.append("#n_families\t").append(candidate.families.length).nl();
				out.append("#encoding\t").append(encoding.name().toLowerCase()).nl();
				out.append("#coordinate_mode\t").append(coordMode.name().toLowerCase()).nl();
				out.append("#min_count_emitted\t").appendA48(minCount).nl();
				out.append("#dense_row_threshold_nonzero_symbols\t").append(DENSE_SYMBOL_THRESHOLD).nl();
				out.append("#start_end_counts\tunavailable").nl();
				out.append("#row_f\tf rank_a48 encoded_rep_id length_a48 consensus_sha80 consensus").nl();
				if(coordMode==CoordMode.EXPLICIT){
					out.append("#row_r\tr anchor_a48 r|d original_sum_a48 symbolCountA48...").nl();
					out.append("#row_i\ti anchor_a48 chain_index_a48 r|d original_sum_a48 symbolCountA48...").nl();
					out.append("#row_d\td anchor_a48 original_del_count_a48").nl();
					out.append("#row_di\tdi anchor_a48 chain_index_a48 r|d original_sum_a48 symbolCountA48...").nl();
				}else{
					out.append("#row_r\tr r|d original_sum_a48 symbolCountA48...").nl();
					out.append("#row_i\ti r|d original_sum_a48 symbolCountA48...").nl();
					out.append("#row_d\td original_del_count_a48").nl();
					out.append("#row_di\tdi r|d original_sum_a48 symbolCountA48...").nl();
				}
				out.append("#symbol_count\tone token: residue letter immediately followed by A48 count").nl();
				for(int i=0; i<candidate.families.length; i++){
					writeFamily(out, candidate.families[i], original.families[i], i,
						encoding, coordMode, minCount);
				}
				out.append("z").nl();
			}
			complete=true;
		}finally{
			if(!complete){Files.deleteIfExists(output);}
		}
	}

	private static void writeFamily(final TextOut out, final Family candidate,
			final Family original, final int rank, final Encoding encoding,
			final CoordMode coordMode, final int minCount) throws IOException{
		out.append('f').tab().appendA48(rank).tab();
		out.appendUtf8(IdentifierCodec.encode(candidate.repId)).tab();
		out.appendA48(candidate.consensus.length).tab();
		out.append(DigestSuffix.fromDigest(candidate.consensusSha256)).tab();
		appendConsensus(out, candidate.consensus);
		out.nl();
		for(int i=0; i<candidate.ref.length; i++){
			writeCounts(out, "r", i, -1, candidate.ref[i], original.ref[i],
				encoding, coordMode, minCount);
			writeChain(out, "i", i, candidate.refIns[i], original.refIns[i],
				encoding, coordMode, minCount);
			if(candidate.del[i].sum>0){
				if(candidate.del[i].sum!=original.del[i].sum){
					throw new IllegalArgumentException("DEL count changed without pruning to zero: "+candidate.repId+" @"+i);
				}
				out.append('d').tab();
				if(coordMode==CoordMode.EXPLICIT){out.appendA48(i).tab();}
				out.appendA48(original.del[i].sum).nl();
			}
			writeChain(out, "di", i, candidate.delIns[i], original.delIns[i],
				encoding, coordMode, minCount);
		}
		out.append('e').nl();
	}

	private static void writeChain(final TextOut out, final String op, final int anchor,
			final Counts[] candidate, final Counts[] original, final Encoding encoding,
			final CoordMode coordMode, final int minCount) throws IOException{
		if(candidate.length>original.length){
			throw new IllegalArgumentException(op+" chain grew at anchor "+anchor+": "+candidate.length+">"+original.length);
		}
		for(int i=0; i<candidate.length; i++){
			writeCounts(out, op, anchor, i, candidate[i], original[i],
				encoding, coordMode, minCount);
		}
	}

	private static void writeCounts(final TextOut out, final String op, final int anchor,
			final int chain, final Counts candidate, final Counts original, final Encoding encoding,
			final CoordMode coordMode, final int minCount) throws IOException{
		candidate.assertSurvivesWithin(original, op, anchor, chain);
		final boolean dense=encoding==Encoding.DENSE ||
			(encoding==Encoding.ADAPTIVE && nonzero(candidate.counts, minCount)>=DENSE_SYMBOL_THRESHOLD);
		out.append(op).tab();
		if(coordMode==CoordMode.EXPLICIT){
			out.appendA48(anchor).tab();
			if(chain>=0){out.appendA48(chain).tab();}
		}
		out.append(dense ? 'd' : 'r').tab().appendA48(original.sum);
		if(dense){
			for(int i=0; i<HbmBundleFormat.NAA; i++){
				final int count=candidate.counts[i];
				out.tab().appendA48(count>=minCount ? count : 0);
			}
		}else{
			for(int i=0; i<HbmBundleFormat.NAA; i++){
				final int count=candidate.counts[i];
				if(count>=minCount){out.tab().append(symbol(i)).appendA48(count);}
			}
		}
		out.nl();
	}

	private static int nonzero(final int[] counts, final int minCount){
		int x=0;
		for(final int count : counts){if(count>=minCount){x++;}}
		return x;
	}

	private static char symbol(final int code){
		if(code==Blosum62.X_CODE){return 'X';}
		if(code>=0 && code<AminoAcid.numberToAcid.length){return (char)AminoAcid.numberToAcid[code];}
		throw new IllegalArgumentException("No residue symbol for encoded slot "+code);
	}

	private static void appendConsensus(final TextOut out, final byte[] consensus) throws IOException{
		for(final byte code : consensus){out.append(symbol(code));}
	}

	private static void fresh(final Path output) throws IOException{
		if(Files.exists(output)){throw new IOException("Output already exists: "+output);}
		final Path parent=output.toAbsolutePath().getParent();
		if(parent!=null && !Files.isDirectory(parent)){throw new IOException("Output parent does not exist: "+parent);}
	}

	private static boolean parseBoolean(final String value, final String key){
		if(value.equalsIgnoreCase("t") || value.equalsIgnoreCase("true")){return true;}
		if(value.equalsIgnoreCase("f") || value.equalsIgnoreCase("false")){return false;}
		throw new IllegalArgumentException(key+" must be t or f");
	}

	private static int parseInt(final String value, final String key){
		try{return Integer.parseInt(value);}
		catch(NumberFormatException nfe){throw new IllegalArgumentException(key+" must be an integer: "+value);}
	}

	enum Encoding{
		SPARSE, DENSE, ADAPTIVE;
		static Encoding parse(final String value){
			for(final Encoding e : values()){if(e.name().equalsIgnoreCase(value)){return e;}}
			throw new IllegalArgumentException("encoding must be sparse, dense, or adaptive: "+value);
		}
	}

	enum CoordMode{
		EXPLICIT, IMPLICIT;
		static CoordMode parse(final String value){
			for(final CoordMode mode : values()){if(mode.name().equalsIgnoreCase(value)){return mode;}}
			throw new IllegalArgumentException("coords must be explicit or implicit: "+value);
		}
	}

	private static final class Bundle{
		final Family[] families;
		final String sourceSha80;

		Bundle(final Family[] families, final String sourceSha80){
			this.families=families;
			this.sourceSha80=sourceSha80;
		}

		static Bundle load(final Path input, final Map<String, byte[]> consensus) throws Exception{
			final byte[] sourceDigest=HbmRareStateProbe.digest(input);
			Path decoded=input, tempDir=null;
			try{
				if(HbmA48Codec.isEncoded(input)){
					tempDir=Files.createTempDirectory("hbm-text-dump-");
					decoded=tempDir.resolve("decoded.mqhb");
					HbmA48Codec.decode(input, decoded);
				}
				return loadNative(decoded, consensus, DigestSuffix.fromDigest(sourceDigest));
			}finally{
				if(tempDir!=null){
					Files.deleteIfExists(decoded);
					Files.deleteIfExists(tempDir);
				}
			}
		}

		private static Bundle loadNative(final Path input, final Map<String, byte[]> consensus,
				final String sourceSha80) throws Exception{
			final long fileLen=Files.size(input);
			try(FileChannel ch=FileChannel.open(input, StandardOpenOption.READ)){
				final byte[] header=read(ch, 0, HbmBundleFormat.FIXED_HEADER_LEN, fileLen);
				if(!HbmBundleFormat.magicMatches(header, HbmBundleFormat.OFF_MAGIC)){
					throw new IllegalArgumentException("bad MQHB magic: "+input);
				}
				if(u16(header, HbmBundleFormat.OFF_FORMAT_VERSION)!=HbmBundleFormat.FORMAT_VERSION){
					throw new IllegalArgumentException("unsupported MQHB version: "+input);
				}
				if(u16(header, HbmBundleFormat.OFF_NAA)!=HbmBundleFormat.NAA ||
						u16(header, HbmBundleFormat.OFF_PAD)!=HbmBundleFormat.EXPECT_PAD){
					throw new IllegalArgumentException("unsupported MQHB NAA/pad: "+input);
				}
				final int n=i32(header, HbmBundleFormat.OFF_N_FAMILIES);
				if(n<1){throw new IllegalArgumentException("n_families<1: "+input);}
				final String[] repIds=new String[n];
				final long[] offsets=new long[n], lengths=new long[n], crcs=new long[n];
				long position=HbmBundleFormat.DIRECTORY_START;
				for(int i=0; i<n; i++){
					final int repLen=i32(read(ch, position, 4, fileLen), 0); position+=4;
					if(repLen<1 || repLen>HbmBundleFormat.CAP_REPID){
						throw new IOException("Invalid representative length at rank "+i);
					}
					repIds[i]=HbmBundleLoader.utf8StrictDecode(read(ch, position, repLen, fileLen)); position+=repLen;
					final byte[] entry=read(ch, position, 20, fileLen); position+=20;
					offsets[i]=i64(entry, 0); lengths[i]=i64(entry, 8); crcs[i]=u32(entry, 16);
				}
				final HashSet<String> seen=new HashSet<String>();
				long expected=position;
				for(int i=0; i<n; i++){
					if(!seen.add(repIds[i])){throw new IllegalArgumentException("Duplicate representative: "+repIds[i]);}
					if(offsets[i]!=expected){throw new IOException("Noncontiguous HBM block at rank "+i);}
					if(lengths[i]<36 || lengths[i]>HbmBundleLoader.CAP_BLOCK_BYTES || (lengths[i]-36)%4!=0){
						throw new IOException("Invalid HBM block length at rank "+i);
					}
					if(offsets[i]<0 || offsets[i]>fileLen || lengths[i]>fileLen-offsets[i]){
						throw new IOException("HBM block out of bounds at rank "+i);
					}
					expected+=lengths[i];
				}
				final byte[] trailer=read(ch, expected, HbmBundleFormat.SHA256_LEN+HbmBundleFormat.MAGIC_LEN, fileLen);
				final byte[] prefix=streamingSha(ch, 0, expected);
				for(int i=0; i<prefix.length; i++){
					if(prefix[i]!=trailer[i]){throw new IOException("MQHB prefix digest mismatch: "+input);}
				}
				if(!HbmBundleFormat.magicMatches(trailer, HbmBundleFormat.SHA256_LEN)){
					throw new IOException("MQHB trailer sentinel mismatch: "+input);
				}
				if(fileLen!=expected+HbmBundleFormat.SHA256_LEN+HbmBundleFormat.MAGIC_LEN){
					throw new IOException("Trailing MQHB bytes: "+input);
				}
				final Family[] families=new Family[n];
				for(int i=0; i<n; i++){
					final byte[] block=read(ch, offsets[i], (int)lengths[i], fileLen);
					if(HbmBundleFormat.crc32(block, 0, block.length)!=crcs[i]){
						throw new IOException("Block CRC mismatch at rank "+i);
					}
					final byte[] cons=consensus.get(repIds[i]);
					if(cons==null){throw new IllegalArgumentException("Consensus missing for "+repIds[i]);}
					families[i]=Family.read(repIds[i], cons, block);
				}
				return new Bundle(families, sourceSha80);
			}
		}

		void assertCompatible(final Bundle original){
			if(families.length!=original.families.length){
				throw new IllegalArgumentException("Candidate/original family counts differ");
			}
			for(int i=0; i<families.length; i++){families[i].assertCompatible(original.families[i], i);}
		}
	}

	private static final class Family{
		final String repId;
		final byte[] consensus, consensusSha256;
		final Counts[] ref, del;
		final Counts[][] refIns, delIns;

		Family(final String repId, final byte[] consensus, final byte[] consensusSha256,
				final Counts[] ref, final Counts[] del, final Counts[][] refIns, final Counts[][] delIns){
			this.repId=repId;
			this.consensus=consensus;
			this.consensusSha256=consensusSha256;
			this.ref=ref;
			this.del=del;
			this.refIns=refIns;
			this.delIns=delIns;
		}

		static Family read(final String repId, final byte[] consensus, final byte[] block){
			final ByteBuffer bb=ByteBuffer.wrap(block);
			final int len=bb.getInt();
			if(len!=consensus.length){throw new IllegalArgumentException("Consensus length mismatch for "+repId);}
			final byte[] cSha=new byte[HbmBundleFormat.SHA256_LEN];
			bb.get(cSha);
			if(!Arrays.equals(cSha, HbmBundleFormat.sha256(consensus, 0, consensus.length))){
				throw new IllegalArgumentException("Consensus SHA mismatch for "+repId);
			}
			final Counts[] ref=new Counts[len], del=new Counts[len];
			final Counts[][] refIns=new Counts[len][], delIns=new Counts[len][];
			for(int i=0; i<len; i++){
				ref[i]=Counts.readHistogram(bb, 1, "REF");
				refIns[i]=readChain(bb);
				del[i]=Counts.readDeletion(bb);
				delIns[i]=readChain(bb);
			}
			if(bb.hasRemaining()){throw new IllegalArgumentException("Trailing block bytes for "+repId);}
			return new Family(repId, consensus.clone(), cSha, ref, del, refIns, delIns);
		}

		private static Counts[] readChain(final ByteBuffer bb){
			final int len=bb.getInt();
			if(len<0 || len>HbmBundleFormat.CAP_CHAIN ||
					(long)len*HbmBundleFormat.INS_RECORD_BYTES>bb.remaining()){
				throw new IllegalArgumentException("Invalid insertion chain length: "+len);
			}
			final Counts[] chain=new Counts[len];
			for(int i=0; i<len; i++){chain[i]=Counts.readHistogram(bb, 1, "INS");}
			return chain;
		}

		void assertCompatible(final Family original, final int rank){
			if(!repId.equals(original.repId)){throw new IllegalArgumentException("rep_id differs at rank "+rank);}
			if(!Arrays.equals(consensus, original.consensus)){
				throw new IllegalArgumentException("Consensus differs at rank "+rank+" ("+repId+")");
			}
			if(!Arrays.equals(consensusSha256, original.consensusSha256)){
				throw new IllegalArgumentException("Consensus digest differs at rank "+rank+" ("+repId+")");
			}
			for(int i=0; i<ref.length; i++){
				if(refIns[i].length>original.refIns[i].length){
					throw new IllegalArgumentException("REF insertion chain grew at "+repId+":"+i);
				}
				if(delIns[i].length>original.delIns[i].length){
					throw new IllegalArgumentException("DEL insertion chain grew at "+repId+":"+i);
				}
			}
		}
	}

	private static final class Counts{
		final int sum;
		final int[] counts;

		Counts(final int sum, final int[] counts){
			this.sum=sum;
			this.counts=counts;
		}

		static Counts readDeletion(final ByteBuffer bb){
			final int sum=bb.getInt();
			if(sum<0){throw new IllegalArgumentException("Negative DEL count");}
			return new Counts(sum, null);
		}

		static Counts readHistogram(final ByteBuffer bb, final int minSum, final String label){
			final int sum=bb.getInt();
			if(sum<minSum){throw new IllegalArgumentException(label+" countSum<"+minSum+": "+sum);}
			final int[] counts=new int[HbmBundleFormat.NAA];
			long total=0;
			for(int i=0; i<counts.length; i++){
				final int count=bb.getInt();
				if(count<0){throw new IllegalArgumentException("Negative "+label+" count");}
				counts[i]=count;
				total+=count;
			}
			if(total!=sum){throw new IllegalArgumentException(label+" sum "+sum+" != counts "+total);}
			return new Counts(sum, counts);
		}

		void assertSurvivesWithin(final Counts original, final String op,
				final int anchor, final int chain){
			if(counts==null || original.counts==null){throw new IllegalArgumentException(op+" is not a residue row");}
			for(int i=0; i<counts.length; i++){
				if(counts[i]!=0 && counts[i]!=original.counts[i]){
					throw new IllegalArgumentException(op+" retained count changed at "+anchor+"/"+chain+
						" symbol "+symbol(i)+": "+counts[i]+" != "+original.counts[i]);
				}
			}
		}
	}

	private static byte[] read(final FileChannel channel, final long position,
			final int length, final long fileLen) throws IOException{
		if(position<0 || length<0 || position>fileLen || length>fileLen-position){
			throw new IOException("Bundle read out of bounds");
		}
		final ByteBuffer buffer=ByteBuffer.allocate(length);
		while(buffer.hasRemaining()){
			final int read=channel.read(buffer, position+buffer.position());
			if(read<=0){throw new IOException("No progress reading bundle");}
		}
		return buffer.array();
	}

	private static byte[] streamingSha(final FileChannel channel, final long from,
			final long to) throws Exception{
		final MessageDigest digest=MessageDigest.getInstance("SHA-256");
		final ByteBuffer buffer=ByteBuffer.allocate(1<<16);
		long position=from;
		while(position<to){
			buffer.clear();
			final long remaining=to-position;
			if(remaining<buffer.capacity()){buffer.limit((int)remaining);}
			final int read=channel.read(buffer, position);
			if(read<=0){throw new IOException("No progress hashing bundle");}
			buffer.flip();
			digest.update(buffer);
			position+=read;
		}
		return digest.digest();
	}

	private static int u16(final byte[] a, final int o){return ((a[o]&0xff)<<8)|(a[o+1]&0xff);}
	private static int i32(final byte[] a, final int o){return ((a[o]&0xff)<<24)|((a[o+1]&0xff)<<16)|((a[o+2]&0xff)<<8)|(a[o+3]&0xff);}
	private static long u32(final byte[] a, final int o){return i32(a, o)&0xffffffffL;}
	private static long i64(final byte[] a, final int o){
		long v=0;
		for(int i=0; i<8; i++){v=(v<<8)|(a[o+i]&0xff);}
		return v;
	}

	private static final class TextOut implements Closeable{
		private final OutputStream out;
		private final ByteBuilder bb=new ByteBuilder(1<<20);

		TextOut(final OutputStream out){this.out=out;}
		TextOut append(final String s){bb.append(s); return flushIfNeeded();}
		TextOut append(final char c){bb.append(c); return flushIfNeeded();}
		TextOut append(final int x){bb.append(x); return flushIfNeeded();}
		TextOut appendA48(final int x){bb.appendA48(x); return flushIfNeeded();}
		TextOut appendUtf8(final String s){
			final byte[] bytes=s.getBytes(java.nio.charset.StandardCharsets.UTF_8);
			bb.append(bytes, 0, bytes.length);
			return flushIfNeeded();
		}
		TextOut tab(){bb.tab(); return this;}
		TextOut nl(){bb.nl(); return flushIfNeeded();}
		private TextOut flushIfNeeded(){
			try{
				if(bb.length()>=1<<20){flush();}
				return this;
			}catch(IOException e){throw new RuntimeException(e);}
		}
		private void flush() throws IOException{
			if(bb.length()>0){
				out.write(bb.array, 0, bb.length());
				bb.clear();
			}
		}
		@Override
		public void close() throws IOException{
			flush();
			if(out instanceof BgzfOutputStream){((BgzfOutputStream)out).writeEOF();}
			out.close();
		}
	}
}
