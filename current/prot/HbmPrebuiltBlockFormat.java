package prot;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.Arrays;
import java.util.Locale;

/**
 * Sealed Pass-5A container for one already-validated MQHB family block. The only writer takes a
 * live {@link HbmBundleBuilder.FamilyInput}, validates its graph, and invokes the canonical block
 * serializer exactly once. Readers independently check the container seal and the complete
 * serialized-block grammar without constructing an {@link AAGraph}.
 *
 * @author Yoimiya
 */
public final class HbmPrebuiltBlockFormat {

	private HbmPrebuiltBlockFormat(){}

	private static final byte[] MAGIC={(byte)'M',(byte)'Q',(byte)'B',(byte)'5'};
	public static final int FORMAT_VERSION=1;
	private static final int TRAILER_LEN=HbmBundleFormat.SHA256_LEN+MAGIC.length;
	private static final int HASH_BUFFER_BYTES=1<<20;

	public static String fileName(final int rank){
		if(rank<0){throw new IllegalArgumentException("PASS5_BLOCK_RANK_RANGE: "+rank);}
		return String.format(Locale.ROOT,"block_%05d.mqb5",rank);
	}

	public static final class Sealed {
		public final int rank, memberCount, blockLength;
		public final long nodeCount, blockCrc32;
		public final String prefixSha256, wholeFileSha256, blockSha256;
		private Sealed(final int rank, final int memberCount, final int blockLength, final long nodeCount,
				final long blockCrc32, final String prefixSha256, final String wholeFileSha256,
				final String blockSha256){
			this.rank=rank; this.memberCount=memberCount; this.blockLength=blockLength; this.nodeCount=nodeCount;
			this.blockCrc32=blockCrc32; this.prefixSha256=prefixSha256; this.wholeFileSha256=wholeFileSha256;
			this.blockSha256=blockSha256;
		}
	}

	/** Opaque verified reference. Only {@link #verify} can construct one. */
	static final class Verified {
		final String path, repId, sourceMqf4PrefixSha256, prefixSha256, wholeFileSha256, blockSha256;
		final byte[] consensus;
		final int rank, memberCount, blockLength;
		final long nodeCount, blockCrc32;
		private Verified(final String path, final String repId, final byte[] consensus, final int rank,
				final int memberCount, final int blockLength, final long nodeCount, final long blockCrc32,
				final String sourceMqf4PrefixSha256, final String prefixSha256,
				final String wholeFileSha256, final String blockSha256){
			this.path=path; this.repId=repId; this.consensus=consensus; this.rank=rank;
			this.memberCount=memberCount; this.blockLength=blockLength; this.nodeCount=nodeCount;
			this.blockCrc32=blockCrc32; this.sourceMqf4PrefixSha256=sourceMqf4PrefixSha256;
			this.prefixSha256=prefixSha256; this.wholeFileSha256=wholeFileSha256;
			this.blockSha256=blockSha256;
		}
	}

	/** Validates and serializes the graph before creating {@code path}; an invalid graph leaves no file. */
	public static Sealed writeNew(final String path, final int rank, final int memberCount,
			final String sourceMqf4PrefixSha256, final HbmBundleBuilder.FamilyInput family) throws IOException{
		if(rank<0){throw new IllegalArgumentException("PASS5_BLOCK_RANK_RANGE: "+rank);}
		if(memberCount<1){throw new IllegalArgumentException("PASS5_BLOCK_MEMBER_COUNT_RANGE: "+memberCount);}
		HbmCanonicalParse.checkLowercaseHexSha256(sourceMqf4PrefixSha256,"source MQF4 prefix sha256");
		if(family==null){throw new IllegalArgumentException("PASS5_BLOCK_NULL_FAMILY");}
		final byte[] rep=HbmBundleBuilder.utf8StrictEncode(family.repId);
		if(rep.length<1 || rep.length>HbmBundleFormat.CAP_REPID){throw new IllegalArgumentException("PASS5_BLOCK_REP_ID_LENGTH: "+rep.length);}
		HbmBundleBuilder.validateGraph(family);
		final byte[] block=HbmBundleBuilder.serializeBlock(family);
		if(block.length<1 || block.length>HbmBundleBuilder.CAP_BLOCK_BYTES){throw new IllegalArgumentException("PASS5_BLOCK_LENGTH_RANGE: "+block.length);}
		final Scan scan=scanBlock(block,family.consensusEnc);
		final long crc=HbmBundleFormat.crc32(block,0,block.length);
		final byte[] blockSha=sha256(block,0,block.length);
		final byte[] sourcePrefix=decodeHex(sourceMqf4PrefixSha256);
		final byte[] consensusSha=sha256(family.consensusEnc,0,family.consensusEnc.length);
		final int headerLength=Math.addExact(132,rep.length);
		final ByteBuffer header=ByteBuffer.allocate(headerLength).order(ByteOrder.BIG_ENDIAN);
		header.put(MAGIC).putShort((short)FORMAT_VERSION).putInt(rank).putShort((short)rep.length).put(rep);
		header.put(sourcePrefix).putInt(family.consensusEnc.length).put(consensusSha).putInt(memberCount);
		header.putLong(scan.nodeCount).putInt(block.length).putInt((int)crc).put(blockSha).flip();
		final MessageDigest prefixDigest=newDigest();
		try(FileChannel ch=FileChannel.open(Paths.get(path),StandardOpenOption.CREATE_NEW,
				StandardOpenOption.WRITE,StandardOpenOption.READ)){
			long pos=0;
			prefixDigest.update(header.array(),0,header.remaining());
			writeFully(ch,pos,header); pos+=headerLength;
			prefixDigest.update(block);
			writeFully(ch,pos,ByteBuffer.wrap(block)); pos+=block.length;
			final byte[] prefix=prefixDigest.digest();
			final ByteBuffer trailer=ByteBuffer.allocate(TRAILER_LEN).put(prefix).put(MAGIC); trailer.flip();
			writeFully(ch,pos,trailer); pos+=TRAILER_LEN; ch.truncate(pos); ch.force(true);
		}
		final String prefixHex=toHex(prefixDigestResult(path));
		//prefixDigestResult hashes only the prefix and independently checks the just-written bytes.
		final String whole=HbmMemberIndexFormat.streamingWholeFileSha256(path);
		return new Sealed(rank,memberCount,block.length,scan.nodeCount,crc,prefixHex,whole,toHex(blockSha));
	}

	/** Fully verifies a sealed container against independently supplied family bindings. */
	static Verified verify(final String path, final int expectedRank, final String expectedRepId,
			final byte[] expectedConsensus, final int expectedMemberCount,
			final String expectedSourceMqf4PrefixSha256, final String expectedPrefixSha256,
			final String expectedWholeFileSha256) throws IOException{
		return verifyInternal(path,expectedRank,expectedRepId,expectedConsensus,expectedMemberCount,
			expectedSourceMqf4PrefixSha256,expectedPrefixSha256,expectedWholeFileSha256,false);
	}

	/** Modern external-provenance verifier. The MQB5 container retains its full internal
	 *  digests; caller-supplied anchors and all failure diagnostics remain sha80. */
	static Verified verifySha80(final String path, final int expectedRank, final String expectedRepId,
			final byte[] expectedConsensus, final int expectedMemberCount,
			final String expectedSourceMqf4PrefixSha80, final String expectedPrefixSha80,
			final String expectedWholeFileSha80) throws IOException{
		return verifyInternal(path,expectedRank,expectedRepId,expectedConsensus,expectedMemberCount,
			expectedSourceMqf4PrefixSha80,expectedPrefixSha80,expectedWholeFileSha80,true);
	}

	private static Verified verifyInternal(final String path, final int expectedRank, final String expectedRepId,
			final byte[] expectedConsensus, final int expectedMemberCount,
			final String expectedSourceMqf4Prefix, final String expectedPrefix,
			final String expectedWholeFile, final boolean sha80) throws IOException{
		if(expectedRank<0 || expectedMemberCount<1){throw new IllegalArgumentException("PASS5_BLOCK_EXPECTED_RANGE");}
		final byte[] expectedRep=HbmBundleBuilder.utf8StrictEncode(expectedRepId);
		if(expectedConsensus==null || expectedConsensus.length<1){throw new IllegalArgumentException("PASS5_BLOCK_EXPECTED_CONSENSUS");}
		if(sha80){
			HbmCanonicalParse.checkLowercaseHexSha80(expectedSourceMqf4Prefix,"expected source MQF4 prefix sha80");
			HbmCanonicalParse.checkLowercaseHexSha80(expectedPrefix,"expected MQB5 prefix sha80");
			HbmCanonicalParse.checkLowercaseHexSha80(expectedWholeFile,"expected MQB5 whole-file sha80");
		}else{
			HbmCanonicalParse.checkLowercaseHexSha256(expectedSourceMqf4Prefix,"expected source MQF4 prefix");
			HbmCanonicalParse.checkLowercaseHexSha256(expectedPrefix,"expected MQB5 prefix");
			HbmCanonicalParse.checkLowercaseHexSha256(expectedWholeFile,"expected MQB5 whole-file hash");
		}
		final HbmMemberIndexFormat.InputIdentity identity=HbmMemberIndexFormat.captureInputIdentity(path);
		if(sha80){
			final String real=DigestSuffix.file(path);
			if(!real.equals(expectedWholeFile)){throw new RuntimeException("PASS5_BLOCK_EXTERNAL_WHOLE_SHA80: real="+
				real+" expected="+expectedWholeFile);}
		}
		final String wholeBefore=HbmMemberIndexFormat.streamingWholeFileSha256(path);
		if(!sha80 && !wholeBefore.equals(expectedWholeFile)){throw new RuntimeException("PASS5_BLOCK_EXTERNAL_WHOLE_HASH: "+wholeBefore);}
		final long fileSize=identity.byteCount;
		if(fileSize<132L+1+TRAILER_LEN){throw new RuntimeException("PASS5_BLOCK_FILE_TOO_SHORT: "+fileSize);}
		final long trailerStart=fileSize-TRAILER_LEN;
		final String prefixBefore=sha256Prefix(path,trailerStart);
		if(sha80){
			final String real=suffix(prefixBefore);
			if(!real.equals(expectedPrefix)){throw new RuntimeException("PASS5_BLOCK_EXTERNAL_PREFIX_SHA80: real="+
				real+" expected="+expectedPrefix);}
		}else if(!prefixBefore.equals(expectedPrefix)){throw new RuntimeException("PASS5_BLOCK_EXTERNAL_PREFIX_HASH: "+prefixBefore);}
		final Parsed parsed=parse(path,fileSize,trailerStart);
		if(parsed.rank!=expectedRank){throw new RuntimeException("PASS5_BLOCK_RANK_MISMATCH: "+parsed.rank+" != "+expectedRank);}
		if(!Arrays.equals(parsed.repId,expectedRep)){throw new RuntimeException("PASS5_BLOCK_REP_ID_MISMATCH: "+expectedRank);}
		if(sha80){
			if(!suffix(parsed.sourcePrefixSha256).equals(expectedSourceMqf4Prefix)){
				throw new RuntimeException("PASS5_BLOCK_SOURCE_PREFIX_SHA80_MISMATCH: "+expectedRank);}
		}else if(!parsed.sourcePrefixSha256.equals(expectedSourceMqf4Prefix)){
			throw new RuntimeException("PASS5_BLOCK_SOURCE_PREFIX_MISMATCH: "+expectedRank);}
		if(parsed.memberCount!=expectedMemberCount){throw new RuntimeException("PASS5_BLOCK_MEMBER_COUNT_MISMATCH: "+parsed.memberCount);}
		if(parsed.consensusLength!=expectedConsensus.length ||
				!Arrays.equals(parsed.consensusSha256,sha256(expectedConsensus,0,expectedConsensus.length))){
			throw new RuntimeException("PASS5_BLOCK_CONSENSUS_MISMATCH: "+expectedRank);
		}
		final Scan scan=scanBlock(parsed.block,expectedConsensus);
		if(scan.nodeCount!=parsed.nodeCount){throw new RuntimeException("PASS5_BLOCK_NODE_COUNT_MISMATCH: "+parsed.nodeCount+" != "+scan.nodeCount);}
		final long crc=HbmBundleFormat.crc32(parsed.block,0,parsed.block.length);
		if(crc!=parsed.blockCrc32){throw new RuntimeException("PASS5_BLOCK_CRC_MISMATCH: "+expectedRank);}
		final String blockSha=toHex(sha256(parsed.block,0,parsed.block.length));
		if(!blockSha.equals(parsed.blockSha256)){throw new RuntimeException("PASS5_BLOCK_SHA_MISMATCH: "+expectedRank);}
		if(!prefixBefore.equals(parsed.storedPrefixSha256)){throw new RuntimeException("PASS5_BLOCK_INTERNAL_PREFIX_MISMATCH: "+expectedRank);}
		HbmMemberIndexFormat.requireUnchangedIdentity(identity,path,"PASS5_BLOCK_INPUT_IDENTITY_CHANGED");
		if(sha80 && !DigestSuffix.file(path).equals(expectedWholeFile)){
			throw new RuntimeException("PASS5_BLOCK_INPUT_SHA80_CHANGED_DURING_VERIFY");
		}
		if(!HbmMemberIndexFormat.streamingWholeFileSha256(path).equals(wholeBefore)){
			throw new RuntimeException("PASS5_BLOCK_INPUT_CHANGED_DURING_VERIFY");
		}
		return new Verified(path,expectedRepId,Arrays.copyOf(expectedConsensus,expectedConsensus.length),expectedRank,
			expectedMemberCount,parsed.block.length,scan.nodeCount,crc,parsed.sourcePrefixSha256,
			prefixBefore,wholeBefore,blockSha);
	}

	/** Re-verifies the immutable source and returns one block copy for immediate assembly. */
	static byte[] loadBlock(final Verified verified) throws IOException{
		final HbmMemberIndexFormat.InputIdentity identity=HbmMemberIndexFormat.captureInputIdentity(verified.path);
		if(!HbmMemberIndexFormat.streamingWholeFileSha256(verified.path).equals(verified.wholeFileSha256)){
			throw new RuntimeException("PASS5_BLOCK_COPY_WHOLE_HASH_BEFORE: rank "+verified.rank);}
		final Parsed parsed=parse(verified.path,identity.byteCount,identity.byteCount-TRAILER_LEN);
		final Scan scan=scanBlock(parsed.block,verified.consensus);
		final String blockSha=toHex(sha256(parsed.block,0,parsed.block.length));
		if(parsed.rank!=verified.rank || !Arrays.equals(parsed.repId,HbmBundleBuilder.utf8StrictEncode(verified.repId)) ||
				!parsed.sourcePrefixSha256.equals(verified.sourceMqf4PrefixSha256) ||
				parsed.consensusLength!=verified.consensus.length ||
				!Arrays.equals(parsed.consensusSha256,sha256(verified.consensus,0,verified.consensus.length)) ||
				parsed.memberCount!=verified.memberCount || parsed.block.length!=verified.blockLength ||
				parsed.nodeCount!=verified.nodeCount || scan.nodeCount!=verified.nodeCount ||
				parsed.blockCrc32!=verified.blockCrc32 || !parsed.blockSha256.equals(verified.blockSha256) ||
				!blockSha.equals(verified.blockSha256) || !parsed.storedPrefixSha256.equals(verified.prefixSha256)){
			throw new RuntimeException("PASS5_BLOCK_METADATA_CHANGED: rank "+verified.rank);}
		HbmMemberIndexFormat.requireUnchangedIdentity(identity,verified.path,"PASS5_BLOCK_COPY_IDENTITY_CHANGED");
		if(!HbmMemberIndexFormat.streamingWholeFileSha256(verified.path).equals(verified.wholeFileSha256)){
			throw new RuntimeException("PASS5_BLOCK_COPY_WHOLE_HASH_AFTER: rank "+verified.rank);}
		return parsed.block;
	}

	private static final class Parsed {
		int rank, consensusLength, memberCount; long nodeCount, blockCrc32;
		byte[] repId, consensusSha256, block; String sourcePrefixSha256, blockSha256, storedPrefixSha256;
	}

	private static Parsed parse(final String path, final long fileSize, final long trailerStart) throws IOException{
		final Parsed p=new Parsed();
		try(FileChannel ch=FileChannel.open(Paths.get(path),StandardOpenOption.READ)){
			long pos=0;
			final ByteBuffer fixed=ByteBuffer.allocate(12).order(ByteOrder.BIG_ENDIAN); readFully(ch,pos,fixed); fixed.flip(); pos+=12;
			final byte[] magic=new byte[4]; fixed.get(magic); if(!Arrays.equals(magic,MAGIC)){throw new RuntimeException("PASS5_BLOCK_MAGIC");}
			if((fixed.getShort()&0xffff)!=FORMAT_VERSION){throw new RuntimeException("PASS5_BLOCK_VERSION");}
			p.rank=fixed.getInt(); final int repLen=fixed.getShort()&0xffff;
			if(repLen<1 || repLen>HbmBundleFormat.CAP_REPID){throw new RuntimeException("PASS5_BLOCK_REP_ID_LENGTH: "+repLen);}
			p.repId=new byte[repLen]; readArray(ch,pos,p.repId); pos+=repLen; HbmCanonicalParse.utf8StrictDecode(p.repId);
			final byte[] source=new byte[32]; readArray(ch,pos,source); pos+=32; p.sourcePrefixSha256=toHex(source);
			final ByteBuffer meta=ByteBuffer.allocate(88).order(ByteOrder.BIG_ENDIAN); readFully(ch,pos,meta); meta.flip(); pos+=88;
			p.consensusLength=meta.getInt(); if(p.consensusLength<1 || p.consensusLength>HbmBundleFormat.CAP_L){throw new RuntimeException("PASS5_BLOCK_CONSENSUS_LENGTH: "+p.consensusLength);}
			p.consensusSha256=new byte[32]; meta.get(p.consensusSha256); p.memberCount=meta.getInt();
			if(p.memberCount<1){throw new RuntimeException("PASS5_BLOCK_MEMBER_COUNT: "+p.memberCount);}
			p.nodeCount=meta.getLong();
			if(p.nodeCount<2 || p.nodeCount>HbmBundleFormat.CAP_TOTAL_NODES){throw new RuntimeException("PASS5_BLOCK_NODE_COUNT: "+p.nodeCount);}
			final int blockLen=meta.getInt(); p.blockCrc32=meta.getInt()&0xffffffffL;
			final byte[] blockSha=new byte[32]; meta.get(blockSha); p.blockSha256=toHex(blockSha);
			if(blockLen<1 || blockLen>HbmBundleBuilder.CAP_BLOCK_BYTES || blockLen>trailerStart-pos){
				throw new RuntimeException("PASS5_BLOCK_LENGTH: "+blockLen);}
			p.block=new byte[blockLen]; readArray(ch,pos,p.block); pos+=blockLen;
			if(pos!=trailerStart){throw new RuntimeException("PASS5_BLOCK_FRAMING: parsed="+pos+" trailer="+trailerStart);}
			final ByteBuffer trailer=ByteBuffer.allocate(TRAILER_LEN); readFully(ch,pos,trailer); trailer.flip();
			final byte[] stored=new byte[32]; trailer.get(stored); p.storedPrefixSha256=toHex(stored);
			final byte[] end=new byte[4]; trailer.get(end); if(!Arrays.equals(end,MAGIC)){throw new RuntimeException("PASS5_BLOCK_TRAILER_MAGIC");}
		}
		return p;
	}

	private static final class Scan {final long nodeCount; Scan(final long nodeCount){this.nodeCount=nodeCount;}}

	/** Bounded structural scan of canonical MQHB block bytes; allocates no graph or per-node object. */
	static Scan scanBlock(final byte[] block, final byte[] consensus){
		final Cursor c=new Cursor(block); final int length=c.i32();
		if(length<1 || length>HbmBundleFormat.CAP_L || consensus==null || consensus.length!=length){
			throw new IllegalArgumentException("PASS5_BLOCK_SCAN_LENGTH: "+length);}
		final byte[] consensusSha=c.bytes(32);
		if(!Arrays.equals(consensusSha,sha256(consensus,0,consensus.length))){throw new IllegalArgumentException("PASS5_BLOCK_SCAN_CONSENSUS_SHA");}
		long nodes=2L*length;
		if(nodes>HbmBundleFormat.CAP_TOTAL_NODES){throw new IllegalArgumentException("PASS5_BLOCK_SCAN_NODE_CAP");}
		for(int pos=0; pos<length; pos++){
			final int refSum=c.i32(); if(refSum<1){throw new IllegalArgumentException("PASS5_BLOCK_SCAN_REF_SUM: "+pos);}
			final long pivotCount=scanCounts(c,refSum,"REF",pos,consensus[pos]&0xff);
			final int residue=consensus[pos]&0xff;
			if(residue>=HbmBundleFormat.NAA || residue==20 || pivotCount<1){throw new IllegalArgumentException("PASS5_BLOCK_SCAN_REF_PIVOT: "+pos);}
			nodes=scanChain(c,nodes,"REF",pos);
			final int delSum=c.i32(); if(delSum<0){throw new IllegalArgumentException("PASS5_BLOCK_SCAN_DEL_SUM: "+pos);}
			nodes=scanChain(c,nodes,"DEL",pos);
		}
		c.end(); return new Scan(nodes);
	}

	private static long scanCounts(final Cursor c, final int expected, final String type, final int pos, final int wantedIndex){
		long sum=0, wanted=-1;
		for(int i=0; i<HbmBundleFormat.NAA; i++){final int v=c.i32();
			if(v<0){throw new IllegalArgumentException("PASS5_BLOCK_SCAN_NEGATIVE_"+type+": "+pos);}
			if(i==wantedIndex){wanted=v;} sum+=v;}
		if(sum!=expected){throw new IllegalArgumentException("PASS5_BLOCK_SCAN_COUNT_SUM_"+type+": "+pos);}
		return wanted;
	}

	private static long scanChain(final Cursor c, long nodes, final String type, final int pos){
		final int length=c.i32();
		if(length<0 || length>HbmBundleFormat.CAP_CHAIN){throw new IllegalArgumentException("PASS5_BLOCK_SCAN_CHAIN_LENGTH_"+type+": "+pos);}
		final long need=(long)length*HbmBundleFormat.INS_RECORD_BYTES;
		if(need>c.remaining()){throw new IllegalArgumentException("PASS5_BLOCK_SCAN_CHAIN_FRAMING_"+type+": "+pos);}
		if(length>HbmBundleFormat.CAP_TOTAL_NODES-nodes){throw new IllegalArgumentException("PASS5_BLOCK_SCAN_NODE_CAP");}
		for(int i=0; i<length; i++){final int sum=c.i32(); if(sum<1){throw new IllegalArgumentException("PASS5_BLOCK_SCAN_INS_SUM");} scanCounts(c,sum,"INS",pos,-1);}
		return nodes+length;
	}

	private static final class Cursor {
		final byte[] a; int p=0; Cursor(final byte[] a){this.a=a;}
		int remaining(){return a.length-p;}
		void need(final int n){if(n<0 || n>a.length-p){throw new IllegalArgumentException("PASS5_BLOCK_SCAN_EOF: need="+n+" pos="+p+" len="+a.length);}}
		int i32(){need(4); final int v=((a[p]&255)<<24)|((a[p+1]&255)<<16)|((a[p+2]&255)<<8)|(a[p+3]&255); p+=4; return v;}
		byte[] bytes(final int n){need(n); final byte[] out=Arrays.copyOfRange(a,p,p+n); p+=n; return out;}
		void end(){if(p!=a.length){throw new IllegalArgumentException("PASS5_BLOCK_SCAN_TRAILING: "+(a.length-p));}}
	}

	private static byte[] prefixDigestResult(final String path) throws IOException{
		final long size=HbmMemberIndexFormat.captureInputIdentity(path).byteCount;
		return decodeHex(sha256Prefix(path,size-TRAILER_LEN));
	}

	private static String sha256Prefix(final String path, final long length) throws IOException{
		final MessageDigest md=newDigest(); final ByteBuffer buffer=ByteBuffer.allocate(HASH_BUFFER_BYTES);
		try(FileChannel ch=FileChannel.open(Paths.get(path),StandardOpenOption.READ)){
			long pos=0;
			while(pos<length){buffer.clear(); buffer.limit((int)Math.min((long)buffer.capacity(),length-pos));
				final int n=ch.read(buffer,pos); if(n<0){throw new IOException("PASS5_BLOCK_HASH_EOF");}
				if(n==0){throw new IOException("PASS5_BLOCK_HASH_NO_PROGRESS");} buffer.flip(); md.update(buffer); pos+=n;}
		}
		return toHex(md.digest());
	}

	private static MessageDigest newDigest(){try{return MessageDigest.getInstance("SHA-256");}catch(NoSuchAlgorithmException e){throw new RuntimeException(e);}}
	private static byte[] sha256(final byte[] a, final int off, final int len){final MessageDigest md=newDigest(); md.update(a,off,len); return md.digest();}
	private static byte[] decodeHex(final String s){HbmCanonicalParse.checkLowercaseHexSha256(s,"sha256"); final byte[] out=new byte[32];
		for(int i=0; i<out.length; i++){out[i]=(byte)((Character.digit(s.charAt(2*i),16)<<4)|Character.digit(s.charAt(2*i+1),16));} return out;}
	private static String suffix(final String full){
		HbmCanonicalParse.checkLowercaseHexSha256(full,"internal sha256");
		return full.substring(full.length()-DigestSuffix.HEX_LENGTH);
	}
	private static String toHex(final byte[] a){return HbmMemberIndexFormat.toHexLower(a);}
	private static void readArray(final FileChannel ch, final long pos, final byte[] a) throws IOException{readFully(ch,pos,ByteBuffer.wrap(a));}
	private static void readFully(final FileChannel ch, final long pos, final ByteBuffer b) throws IOException{long p=pos; while(b.hasRemaining()){
		final int n=ch.read(b,p); if(n<0){throw new IOException("PASS5_BLOCK_EOF: "+p);} if(n==0){throw new IOException("PASS5_BLOCK_READ_NO_PROGRESS");} p+=n;}}
	private static void writeFully(final FileChannel ch, final long pos, final ByteBuffer b) throws IOException{long p=pos; while(b.hasRemaining()){
		final int n=ch.write(b,p); if(n<=0){throw new IOException("PASS5_BLOCK_WRITE_NO_PROGRESS");} p+=n;}}
}
