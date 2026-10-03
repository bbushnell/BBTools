package prot;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.CharBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.CharacterCodingException;
import java.nio.charset.CharsetEncoder;
import java.nio.charset.CodingErrorAction;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.HashSet;
import java.util.List;

import structures.ByteBuilder;

/**
 * Streaming builder for the deterministic HBM model-bundle (design v5). ALL families are fully
 * validated (canonical knobs raw-bit-exact; every REF/DEL/INS node field; {@code weight==count} and
 * {@code weightSum==countSum}; REF pivot count &ge; 1; valid consensus residues; caps) BEFORE the
 * output file is created, so a malformed later family never leaves a partial artifact. Writes to a
 * {@link FileChannel} with a reserved directory back-filled per family, holding at most one block in
 * memory and never casting a whole-file length to {@code int} (offsets/lengths are {@code long}).
 *
 * <p><b>Publication is external.</b> Java-8 NIO cannot guarantee an atomic <i>no-replace</i>
 * directory rename (target-exists behavior of {@code ATOMIC_MOVE} is provider-specific), so this
 * class does NOT publish. The Bash orchestrator stages this bundle in a fresh temp directory and
 * publishes it with its tested no-replace/state-verification pattern.</p>
 *
 * <p>Scope: serialization only — no integration, no real corpus run, no HBM form/cutoff/default.</p>
 *
 * @author UMP45
 */
public final class HbmBundleBuilder {

	private HbmBundleBuilder(){}

	/** Writes a lossless A48/BGZF representation of an existing native bundle to a fresh path. */
	public static void encodeA48(final Path nativeBundle, final Path output) throws IOException{
		HbmA48Codec.encode(nativeBundle, output);
	}

	interface PrebuiltCopyHook {void beforeCopy(int rank, String path) throws IOException;}
	static PrebuiltCopyHook prebuiltCopyHookForTest=null;

	/** Documented per-block serialized cap (1 GiB); must match {@code HbmBundleLoader.CAP_BLOCK_BYTES}. */
	public static final int CAP_BLOCK_BYTES=1<<30;

	public static final class FamilyInput{
		public final String repId; public final byte[] consensusEnc; public final AAGraph graph;
		public FamilyInput(String repId, byte[] consensusEnc, AAGraph graph){
			if(repId==null||repId.isEmpty()){throw new IllegalArgumentException("empty repId");}
			if(consensusEnc==null||consensusEnc.length==0){throw new IllegalArgumentException("empty consensus");}
			if(graph==null){throw new IllegalArgumentException("null graph");}
			this.repId=repId; this.consensusEnc=consensusEnc; this.graph=graph;
		}
	}

	public static void build(final Path out, final List<FamilyInput> families, final byte[][] provenance) throws IOException{
		if(families==null || families.isEmpty()){throw new IllegalArgumentException("nonempty family list required");}
		checkProvenance(provenance);
		final int n=families.size();

		//---- Preflight: validate EVERY family + all caps BEFORE creating the file ----
		final byte[][] repBytes=new byte[n][];
		final long[] blockLens=new long[n];
		final HashSet<String> repSeen=new HashSet<String>();
		long dirLen=0, totalNodes=0, sumBlockBytes=0;
		for(int i=0; i<n; i++){
			final FamilyInput fi=families.get(i);
			if(fi==null){throw new IllegalArgumentException("null FamilyInput at "+i);}
			if(!repSeen.add(fi.repId)){throw new IllegalArgumentException("duplicate rep_id: "+fi.repId);}
			repBytes[i]=utf8StrictEncode(fi.repId);
			if(repBytes[i].length<1 || repBytes[i].length>HbmBundleFormat.CAP_REPID){throw new IllegalArgumentException("rep_id_len out of range at "+i);}
			validateGraph(fi);
			final long bl=computeBlockLen(fi);
			if(bl>CAP_BLOCK_BYTES){throw new IllegalArgumentException("block "+i+" exceeds CAP_BLOCK_BYTES");}
			blockLens[i]=bl;
			final long nodes=2L*fi.graph.ref.length+chainNodeCount(fi);
			if(nodes>HbmBundleFormat.CAP_TOTAL_NODES-totalNodes){throw new IllegalArgumentException("CAP_TOTAL_NODES exceeded");}
			totalNodes+=nodes;
			if(bl>HbmBundleFormat.CAP_TOTAL_BYTES-sumBlockBytes){throw new IllegalArgumentException("CAP_TOTAL_BYTES exceeded");}
			sumBlockBytes+=bl;
			dirLen+=HbmBundleFormat.DIR_ENTRY_FIXED_BYTES+repBytes[i].length;
		}
		final long directoryStart=HbmBundleFormat.DIRECTORY_START;
		final long blocksStart=directoryStart+dirLen;

		try(FileChannel ch=FileChannel.open(out, StandardOpenOption.CREATE_NEW, StandardOpenOption.WRITE, StandardOpenOption.READ)){
			writeHeader(ch, n, provenance);
			for(int i=0; i<n; i++){//directory: rep_id + 20-byte zero placeholder
				final ByteBuilder e=new ByteBuilder();
				putI32(e, repBytes[i].length); e.append(repBytes[i], 0, repBytes[i].length);
				putI64(e, 0); putI64(e, 0); putU32(e, 0);
				writeFully(ch, e);
			}
			long off=blocksStart;
			for(int i=0; i<n; i++){
				final byte[] block=serializeBlock(families.get(i));
				if(block.length!=blockLens[i]){throw new IllegalStateException("block "+i+" length "+block.length+" != preflight "+blockLens[i]);}
				final long crc=HbmBundleFormat.crc32(block, 0, block.length);
				ch.position(off); writeFully(ch, ByteBuffer.wrap(block));
				backfillDirEntry(ch, repBytes, i, off, block.length, crc);
				off+=block.length;
			}
			final long trailerStart=off;
			final byte[] prefixSha=streamingSha(ch, 0, trailerStart);
			ch.position(trailerStart);
			final ByteBuilder t=new ByteBuilder();
			t.append(prefixSha, 0, prefixSha.length);
			t.append(HbmBundleFormat.magic(), 0, HbmBundleFormat.MAGIC_LEN);
			writeFully(ch, t);
			ch.truncate(trailerStart+HbmBundleFormat.SHA256_LEN+HbmBundleFormat.MAGIC_LEN);
		}
	}

	/**
	 * Assembles the existing MQHB byte format from sealed, prevalidated family blocks without
	 * reconstructing a graph or serializing a block. The opaque references can only be produced by
	 * {@link HbmPrebuiltBlockFormat#verify}; every source is reverified before output creation and
	 * again while its exact block bytes are copied, closing the preflight/copy mutation window.
	 */
	static void assemblePrebuilt(final Path out, final List<HbmPrebuiltBlockFormat.Verified> supplied,
			final byte[][] provenance) throws IOException{
		if(supplied==null || supplied.isEmpty()){throw new IllegalArgumentException("nonempty prebuilt block list required");}
		checkProvenance(provenance);
		final int n=supplied.size();
		final HbmPrebuiltBlockFormat.Verified[] blocks=new HbmPrebuiltBlockFormat.Verified[n];
		final byte[][] repBytes=new byte[n][];
		final HashSet<String> repSeen=new HashSet<String>();
		long dirLen=0, totalNodes=0, sumBlockBytes=0;
		for(int i=0; i<n; i++){
			final HbmPrebuiltBlockFormat.Verified v=supplied.get(i);
			if(v==null){throw new IllegalArgumentException("null prebuilt block at "+i);}
			if(v.rank!=i){throw new IllegalArgumentException("prebuilt rank/order mismatch: "+v.rank+" != "+i);}
			final HbmPrebuiltBlockFormat.Verified fresh=HbmPrebuiltBlockFormat.verify(v.path,v.rank,v.repId,
				v.consensus,v.memberCount,v.sourceMqf4PrefixSha256,v.prefixSha256,v.wholeFileSha256);
			blocks[i]=fresh;
			if(!repSeen.add(fresh.repId)){throw new IllegalArgumentException("duplicate rep_id: "+fresh.repId);}
			repBytes[i]=utf8StrictEncode(fresh.repId);
			if(repBytes[i].length<1 || repBytes[i].length>HbmBundleFormat.CAP_REPID){throw new IllegalArgumentException("rep_id_len out of range at "+i);}
			if(fresh.blockLength<1 || fresh.blockLength>CAP_BLOCK_BYTES){throw new IllegalArgumentException("block "+i+" exceeds CAP_BLOCK_BYTES");}
			if(fresh.nodeCount>HbmBundleFormat.CAP_TOTAL_NODES-totalNodes){throw new IllegalArgumentException("CAP_TOTAL_NODES exceeded");}
			totalNodes=Math.addExact(totalNodes,fresh.nodeCount);
			if(fresh.blockLength>HbmBundleFormat.CAP_TOTAL_BYTES-sumBlockBytes){throw new IllegalArgumentException("CAP_TOTAL_BYTES exceeded");}
			sumBlockBytes=Math.addExact(sumBlockBytes,fresh.blockLength);
			dirLen=Math.addExact(dirLen,HbmBundleFormat.DIR_ENTRY_FIXED_BYTES+repBytes[i].length);
		}
		final long blocksStart=Math.addExact((long)HbmBundleFormat.DIRECTORY_START,dirLen);

		try(FileChannel ch=FileChannel.open(out,StandardOpenOption.CREATE_NEW,
				StandardOpenOption.WRITE,StandardOpenOption.READ)){
			writeHeader(ch,n,provenance);
			for(int i=0; i<n; i++){
				final ByteBuilder entry=new ByteBuilder();
				putI32(entry,repBytes[i].length); entry.append(repBytes[i],0,repBytes[i].length);
				putI64(entry,0); putI64(entry,0); putU32(entry,0); writeFully(ch,entry);
			}
			long off=blocksStart;
			for(int i=0; i<n; i++){
				final PrebuiltCopyHook hook=prebuiltCopyHookForTest;
				if(hook!=null){hook.beforeCopy(i,blocks[i].path);}
				final byte[] block=HbmPrebuiltBlockFormat.loadBlock(blocks[i]);
				if(block.length!=blocks[i].blockLength){throw new IllegalStateException("prebuilt block length changed at "+i);}
				final long crc=HbmBundleFormat.crc32(block,0,block.length);
				if(crc!=blocks[i].blockCrc32){throw new IllegalStateException("prebuilt block CRC changed at "+i);}
				ch.position(off); writeFully(ch,ByteBuffer.wrap(block));
				backfillDirEntry(ch,repBytes,i,off,block.length,crc);
				off=Math.addExact(off,block.length);
			}
			final long trailerStart=off;
			final byte[] prefixSha=streamingSha(ch,0,trailerStart);
			ch.position(trailerStart);
			final ByteBuilder trailer=new ByteBuilder();
			trailer.append(prefixSha,0,prefixSha.length);
			trailer.append(HbmBundleFormat.magic(),0,HbmBundleFormat.MAGIC_LEN);
			writeFully(ch,trailer);
			ch.truncate(Math.addExact(trailerStart,HbmBundleFormat.SHA256_LEN+HbmBundleFormat.MAGIC_LEN));
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------      Graph validation         ----------------*/
	/*--------------------------------------------------------------*/

	/** Proves the family's graph is a canonical, self-consistent scoring model before any byte is written. */
	static void validateGraph(final FamilyInput fi){
		final AAGraph g=fi.graph;
		if(g.pad!=HbmBundleFormat.EXPECT_PAD){throw new IllegalArgumentException("graph pad!=0");}
		if(g.pivot.length!=fi.consensusEnc.length){throw new IllegalArgumentException("pivot length != consensus");}
		final int L=g.ref.length;
		if(L<1 || L>HbmBundleFormat.CAP_L){throw new IllegalArgumentException("L out of range: "+L);}
		if(g.del.length!=L){throw new IllegalArgumentException("del length != ref length");}
		for(int i=0; i<L; i++){
			if(g.pivot[i]!=fi.consensusEnc[i]){throw new IllegalArgumentException("pivot != consensus at "+i);}
			if(!isValidResidue(fi.consensusEnc[i])){throw new IllegalArgumentException("invalid consensus residue at "+i+": "+fi.consensusEnc[i]);}
		}
		//Knobs raw-bit-exact.
		if(g.weightByIdentity!=HbmBundleFormat.KNOB_WEIGHT_BY_IDENTITY || g.baseWeight!=HbmBundleFormat.KNOB_BASE_WEIGHT
				|| g.minDepth!=HbmBundleFormat.KNOB_MIN_DEPTH
				|| fb(g.MAF_sub)!=fb(HbmBundleFormat.KNOB_MAF_SUB) || fb(g.MAF_del)!=fb(HbmBundleFormat.KNOB_MAF_DEL)
				|| fb(g.MAF_ins)!=fb(HbmBundleFormat.KNOB_MAF_INS) || fb(g.identityCeiling)!=fb(HbmBundleFormat.KNOB_IDENTITY_CEILING)
				|| fb(g.trimDepthFraction)!=fb(HbmBundleFormat.KNOB_TRIM_DEPTH_FRACTION)){
			throw new IllegalArgumentException("graph knobs are not canonical");
		}
		for(int i=0; i<L; i++){
			final AAGraphNode ref=g.ref[i];
			if(ref.type!=AAGraphNode.REF || ref.rpos!=i || ref.refResidue!=g.pivot[i]){throw new IllegalArgumentException("REF node meta @"+i);}
			validateHist(ref, false);
			if(ref.countSum<1){throw new IllegalArgumentException("REF countSum<1 @"+i);}
			if(ref.count[g.pivot[i]&0xff]<1){throw new IllegalArgumentException("REF pivot count<1 @"+i);}
			validateChain(ref.insEdge, i+1);
			final AAGraphNode del=g.del[i];
			if(del.type!=AAGraphNode.DEL || del.rpos!=i || del.refResidue!=g.pivot[i]){throw new IllegalArgumentException("DEL node meta @"+i);}
			if(del.count!=null || del.weight!=null){throw new IllegalArgumentException("DEL arrays must be null @"+i);}
			if(del.countSum<0 || del.weightSum!=del.countSum){throw new IllegalArgumentException("DEL sums @"+i);}
			validateChain(del.insEdge, i+1);
		}
	}
	private static void validateChain(AAGraphNode ins, final int rpos){
		int len=0;
		for(AAGraphNode p=ins; p!=null; p=p.insEdge){
			len++; if(len>HbmBundleFormat.CAP_CHAIN){throw new IllegalArgumentException("chain exceeds CAP_CHAIN");}
			if(p.type!=AAGraphNode.INS || p.rpos!=rpos || p.refResidue!=Blosum62.X_CODE){throw new IllegalArgumentException("INS node meta");}
			validateHist(p, false);
			if(p.countSum<1){throw new IllegalArgumentException("INS countSum<1");}
		}
	}
	/** Checks a REF/INS node's histograms: length NAA, non-negative, weight==count, weightSum==countSum==Σcount. */
	private static void validateHist(final AAGraphNode node, final boolean del){
		if(node.count==null || node.weight==null || node.count.length!=HbmBundleFormat.NAA || node.weight.length!=HbmBundleFormat.NAA){
			throw new IllegalArgumentException("node histogram length/nullity");
		}
		long sum=0;
		for(int k=0; k<HbmBundleFormat.NAA; k++){
			if(node.count[k]<0 || node.weight[k]<0){throw new IllegalArgumentException("negative count/weight");}
			if(node.weight[k]!=node.count[k]){throw new IllegalArgumentException("weight[k]!=count[k] (non-canonical build)");}
			sum+=node.count[k];
		}
		if(sum!=node.countSum || node.weightSum!=node.countSum){throw new IllegalArgumentException("countSum/weightSum != Σcount");}
	}
	private static int fb(final float f){return Float.floatToRawIntBits(f);}

	/** True iff {@code b} is a valid encoded protein residue for a consensus/pivot byte: real residues
	 *  0-19, or the unknown code {@link Blosum62#X_CODE} (21). Encoded stop (20 — deliberately unused
	 *  by {@link AAGraphNode}'s histogram slots, rejected by {@code Blosum62.encode}) and any other
	 *  value (negative, or above 21) are invalid. Duplicated identically in {@code HbmBundleLoader}
	 *  rather than promoted into the already-accepted-and-applied {@code HbmBundleFormat} (Increment 1)
	 *  — see the v3→v4 delta doc if this should move there instead. */
	private static boolean isValidResidue(final byte b){
		return (b>=0 && b<=19) || b==Blosum62.X_CODE;
	}

	private static long computeBlockLen(final FamilyInput fi){
		final AAGraph g=fi.graph; final int L=g.ref.length;
		long len=4+HbmBundleFormat.SHA256_LEN;//consensus_len + consensus_sha256
		for(int i=0; i<L; i++){
			len+=4+HbmBundleFormat.COUNT_ARRAY_BYTES;//ref countSum + count[22]
			len+=4+(long)chainLen(g.ref[i].insEdge)*HbmBundleFormat.INS_RECORD_BYTES;//ref chain_len + nodes
			len+=4;//del countSum
			len+=4+(long)chainLen(g.del[i].insEdge)*HbmBundleFormat.INS_RECORD_BYTES;//del chain_len + nodes
		}
		return len;
	}
	private static long chainNodeCount(final FamilyInput fi){
		final AAGraph g=fi.graph; long c=0;
		for(int i=0; i<g.ref.length; i++){c+=chainLen(g.ref[i].insEdge)+chainLen(g.del[i].insEdge);}
		return c;
	}
	private static int chainLen(AAGraphNode n){int c=0; for(; n!=null; n=n.insEdge){c++;} return c;}

	/*--------------------------------------------------------------*/
	/*----------------        Serialization          ----------------*/
	/*--------------------------------------------------------------*/

	static byte[] serializeBlock(final FamilyInput fi){
		final AAGraph g=fi.graph; final int L=g.ref.length;
		final ByteBuilder b=new ByteBuilder();
		putI32(b, L);
		final byte[] cSha=HbmBundleFormat.sha256(fi.consensusEnc, 0, fi.consensusEnc.length);
		b.append(cSha, 0, cSha.length);
		for(int i=0; i<L; i++){
			putI32(b, g.ref[i].countSum); writeCount(b, g.ref[i]); writeInsChain(b, g.ref[i].insEdge);
			putI32(b, g.del[i].countSum); writeInsChain(b, g.del[i].insEdge);
		}
		return b.toBytes();
	}
	private static void writeCount(final ByteBuilder b, final AAGraphNode node){for(int k=0; k<HbmBundleFormat.NAA; k++){putI32(b, node.count[k]);}}
	private static void writeInsChain(final ByteBuilder b, AAGraphNode ins){
		putI32(b, chainLen(ins));
		for(AAGraphNode p=ins; p!=null; p=p.insEdge){putI32(b, p.countSum); writeCount(b, p);}
	}

	/*--------------------------------------------------------------*/
	/*----------------      Low-level writers        ----------------*/
	/*--------------------------------------------------------------*/

	private static void writeHeader(final FileChannel ch, final int n, final byte[][] provenance) throws IOException{
		final ByteBuilder h=new ByteBuilder();
		h.append(HbmBundleFormat.magic(), 0, HbmBundleFormat.MAGIC_LEN);
		putU16(h, HbmBundleFormat.FORMAT_VERSION); putU16(h, HbmBundleFormat.NAA); putU16(h, HbmBundleFormat.EXPECT_PAD); putI32(h, n);
		h.append((byte)(HbmBundleFormat.KNOB_WEIGHT_BY_IDENTITY?1:0));
		putI32(h, HbmBundleFormat.KNOB_BASE_WEIGHT);
		putI32(h, fb(HbmBundleFormat.KNOB_MAF_SUB)); putI32(h, fb(HbmBundleFormat.KNOB_MAF_DEL)); putI32(h, fb(HbmBundleFormat.KNOB_MAF_INS));
		putI32(h, HbmBundleFormat.KNOB_MIN_DEPTH); putI32(h, fb(HbmBundleFormat.KNOB_IDENTITY_CEILING)); putI32(h, fb(HbmBundleFormat.KNOB_TRIM_DEPTH_FRACTION));
		for(int i=0; i<HbmBundleFormat.PROVENANCE_COUNT; i++){h.append(provenance[i], 0, HbmBundleFormat.SHA256_LEN);}
		if(h.length()!=HbmBundleFormat.FIXED_HEADER_LEN){throw new IllegalStateException("header len "+h.length()+" != "+HbmBundleFormat.FIXED_HEADER_LEN);}
		ch.position(0); writeFully(ch, h);
	}
	private static void checkProvenance(final byte[][] provenance){
		if(provenance==null || provenance.length!=HbmBundleFormat.PROVENANCE_COUNT){throw new IllegalArgumentException("provenance must have "+HbmBundleFormat.PROVENANCE_COUNT+" fields");}
		for(int i=0; i<provenance.length; i++){if(provenance[i]==null||provenance[i].length!=HbmBundleFormat.SHA256_LEN){throw new IllegalArgumentException("provenance["+i+"] must be 32 bytes");}}
		if(allZero(provenance[HbmBundleFormat.PROV_LOADER_SOURCE]) || allZero(provenance[HbmBundleFormat.PROV_LOADER_CLASS])){
			throw new IllegalArgumentException("loader_source/class provenance must be nonzero in a built bundle");
		}
	}
	private static boolean allZero(final byte[] a){for(byte b : a){if(b!=0){return false;}} return true;}

	private static long dirEntryStart(final byte[][] repBytes, final int idx){
		long p=HbmBundleFormat.DIRECTORY_START; for(int i=0; i<idx; i++){p+=HbmBundleFormat.DIR_ENTRY_FIXED_BYTES+repBytes[i].length;} return p;
	}
	private static void backfillDirEntry(final FileChannel ch, final byte[][] repBytes, final int idx,
			final long blockOffset, final long blockLen, final long crc) throws IOException{
		final long triplePos=dirEntryStart(repBytes, idx)+4+repBytes[idx].length;
		final ByteBuffer bb=ByteBuffer.allocate(HbmBundleFormat.DIR_ENTRY_FIXED_BYTES-4);
		bb.putLong(blockOffset); bb.putLong(blockLen); bb.putInt((int)(crc & 0xffffffffL)); bb.flip();
		ch.position(triplePos); while(bb.hasRemaining()){ch.write(bb);}
	}
	private static byte[] streamingSha(final FileChannel ch, final long from, final long to) throws IOException{
		final MessageDigest md; try{md=MessageDigest.getInstance("SHA-256");}catch(NoSuchAlgorithmException e){throw new RuntimeException("SHA-256 unavailable", e);}
		final ByteBuffer buf=ByteBuffer.allocate(1<<16); long pos=from;
		while(pos<to){buf.clear(); final long rem=to-pos; if(rem<buf.capacity()){buf.limit((int)rem);}
			final int r=ch.read(buf, pos); if(r<0){throw new IOException("EOF hashing prefix at "+pos);} buf.flip(); md.update(buf); pos+=r;}
		return md.digest();
	}
	private static void putU16(final ByteBuilder b, final int v){b.append((byte)((v>>>8)&0xff)); b.append((byte)(v&0xff));}
	private static void putI32(final ByteBuilder b, final int v){b.append((byte)((v>>>24)&0xff)); b.append((byte)((v>>>16)&0xff)); b.append((byte)((v>>>8)&0xff)); b.append((byte)(v&0xff));}
	private static void putU32(final ByteBuilder b, final long v){putI32(b, (int)(v & 0xffffffffL));}
	private static void putI64(final ByteBuilder b, final long v){for(int s=56; s>=0; s-=8){b.append((byte)((v>>>s)&0xff));}}
	private static void writeFully(final FileChannel ch, final ByteBuilder b) throws IOException{writeFully(ch, ByteBuffer.wrap(b.toBytes()));}
	private static void writeFully(final FileChannel ch, final ByteBuffer bb) throws IOException{while(bb.hasRemaining()){ch.write(bb);}}

	static byte[] utf8StrictEncode(final String s){
		final CharsetEncoder enc=StandardCharsets.UTF_8.newEncoder().onMalformedInput(CodingErrorAction.REPORT).onUnmappableCharacter(CodingErrorAction.REPORT);
		try{final ByteBuffer bb=enc.encode(CharBuffer.wrap(s)); final byte[] out=new byte[bb.remaining()]; bb.get(out); return out;}
		catch(CharacterCodingException e){throw new IllegalArgumentException("rep_id not valid UTF-16/UTF-8: "+e.getMessage());}
	}
}
