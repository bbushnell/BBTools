package prot;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.CharBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.CharacterCodingException;
import java.nio.charset.CharsetDecoder;
import java.nio.charset.CharsetEncoder;
import java.nio.charset.CodingErrorAction;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import fileIO.ByteFile;

/**
 * Eager loader for the deterministic HBM model-bundle (design v5). Reads over a {@link FileChannel}
 * with {@code long} positions (no whole-file {@code int} cast) and reconstructs each family's
 * {@link AAGraph} into a PRIVATE, frozen model. Every read is remaining-byte checked; the
 * runtime-semantic hash check and all-zero-loader-hash rejection run BEFORE any graph is built; the
 * total-node and byte caps, and consensus-residue validity, are enforced BEFORE node allocation.
 *
 * <p>{@link Loaded} exposes no graph/node/array/pivot alias — only family count, immutable rep IDs,
 * and a scoring method. Structural inspection for tests is a package-private compare-against-expected
 * that never returns the private graph.</p>
 *
 * <p><b>v3:</b> {@link Loaded#score} forwards {@code passedAllOtherFilters} instead of hardcoding
 * {@code true}.</p>
 *
 * <p><b>v4 (root_review_hbm_increment2_v3, sha256 49c83ea0...):</b> (1) a cloned provider consensus
 * containing an invalid residue (encoded stop 20, negative, or above {@link Blosum62#X_CODE}) is now
 * rejected with an explicit message BEFORE {@link AAGraph} construction, instead of the array index
 * either silently succeeding (stop 20, in-bounds) or crashing with an unexplained
 * {@link ArrayIndexOutOfBoundsException} (out-of-range bytes) from inside {@code AAGraph}'s own
 * scaffold; (2) {@link #readChain} now checks the declared {@code chain_len} against the COMPLETE
 * remaining block bytes before allocating any {@link AAGraphNode}, not just against {@code CAP_CHAIN};
 * (3) the pre-directory-allocation {@code n_families} floor now also reserves one minimum {@code L=1}
 * family block per family plus the 36-byte trailer, not just minimum directory bytes.</p>
 *
 * @author UMP45
 */
public final class HbmBundleLoader {

	private HbmBundleLoader(){}

	/** Documented per-block serialized cap (1 GiB): a single family block never exceeds this, so it
	 *  fits a {@code byte[]} and the cap — not array-allocation failure — is the limit. Must match
	 *  {@code HbmBundleBuilder.CAP_BLOCK_BYTES}. */
	public static final int CAP_BLOCK_BYTES=1<<30;
	private static final String[] SEMANTIC_PROVENANCE_KEYS={
		"format_spec", "loader_source", "loader_class", "aagraph",
		"aagraphnode", "aagraphscorer", "glocalaminolinear", "blosum62",
		"builder_source", "builder_class", "roster", "consensus_ref",
		"source_corpus", "cluster_membership", "member_policy", "member_manifest"
	};

	/** Strictly loads the canonical ordered 16-row runtime/build semantic-provenance table. */
	public static byte[][] loadSemanticProvenance(final String path){
		final byte[][] out=new byte[HbmBundleFormat.PROVENANCE_COUNT][];
		int row=0;
		final ByteFile bf=ByteFile.makeByteFile(path,true);
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(row>=SEMANTIC_PROVENANCE_KEYS.length){throw new IllegalArgumentException("PASS5B_PROVENANCE_EXTRA_ROW: "+row);}
				final String s=HbmCanonicalParse.utf8StrictDecode(line);
				final int tab=s.indexOf('\t'),tab2=(tab<0 ? -1 : s.indexOf('\t',tab+1));
				if(tab<1 || (tab2>=0 && (tab2+1>=s.length() || s.indexOf('\t',tab2+1)>=0))){
					throw new IllegalArgumentException("PASS5B_PROVENANCE_ROW_SHAPE: "+row);
				}
				final String base=SEMANTIC_PROVENANCE_KEYS[row], observedKey=s.substring(0,tab);
				final boolean modern=observedKey.equals(base+"_sha80"), legacy=observedKey.equals(base+"_sha256");
				if(!(modern || legacy)){throw new IllegalArgumentException("PASS5B_PROVENANCE_KEY: "+row);}
				//The accepted production manifest carries an optional third provenance-path
				//column.  The path is explanatory metadata; the trusted semantic value is the
				//middle digest.  Two-column synthetic/legacy manifests remain valid.
				final String text=s.substring(tab+1,tab2<0 ? s.length() : tab2);
				if((modern && text.length()!=DigestSuffix.HEX_LENGTH) || (legacy && text.length()!=64)){
					throw new IllegalArgumentException("PASS5B_PROVENANCE_DIGEST_WIDTH: "+row);
				}
				out[row]=DigestSuffix.decodeSuffix(DigestSuffix.normalizeRecorded(text,observedKey),observedKey);
				row++;
			}
		}finally{bf.close();}
		if(row!=SEMANTIC_PROVENANCE_KEYS.length){throw new IllegalArgumentException("PASS5B_PROVENANCE_ROW_COUNT: "+row);}
		return out;
	}

	/** Minimum possible serialized family block: L=1, zero insertion chains — 4 (L) + SHA256_LEN
	 *  (consensus_sha256) + one column's floor: 4 (ref countSum) + COUNT_ARRAY_BYTES (ref count[]) +
	 *  4 (ref chain_len=0) + 4 (del countSum) + 4 (del chain_len=0). Used only to size the
	 *  pre-directory-allocation {@code n_families} floor below — never assumed elsewhere. */
	private static final long MIN_BLOCK_BYTES=4+HbmBundleFormat.SHA256_LEN+4+HbmBundleFormat.COUNT_ARRAY_BYTES+4+4+4;

	public interface ConsensusProvider{ byte[] consensusFor(String repId); }

	/** The loaded, verified bundle. Graphs are PRIVATE and frozen; nothing here returns an alias to them. */
	public static final class Loaded{
		private final String[] repIds;      //Strings are immutable
		private final AAGraph[] graphs;     //never exposed
		Loaded(String[] repIds, AAGraph[] graphs){this.repIds=repIds; this.graphs=graphs;}
		public int familyCount(){return graphs.length;}
		/** Declares the stable ordering and meaning of values returned by {@link #score}. */
		public String scoreContract(){return HbmScoreContract.VALUE;}
		public String repId(final int i){return repIds[i];}
		/** Derives detached immutable experimental position scores; loaded graphs stay private and unchanged. */
		HbmPositionModel[] positionModels(String kind,double beta,boolean clip){
			return HbmPositionModel.derive(graphs,kind,beta,clip);
		}
		/** As above, with an explicit experimental lower clipping bound in half-bit units. */
		HbmPositionModel[] positionModels(String kind,double beta,boolean clip,double clipMin){
			return HbmPositionModel.derive(graphs,kind,beta,clip,clipMin);
		}
		/** Returns a new background vector; no private histogram storage escapes. */
		double[] positionBackground(){return HbmPositionModel.backgroundProbabilities(graphs);}
		/** Refines one family into fresh caller-owned output while leaving the loaded graph untouched. */
		HbmProfileRefiner.Result refineFamily(int rank,byte[][] members,double[] background,int padding){
			if(rank<0 || rank>=graphs.length){throw new IllegalArgumentException("Refinement rank outside loaded roster: "+rank);}
			return HbmProfileRefiner.refine(graphs[rank],members,background,padding);
		}
		/** Composes verified immutable singleton bundles in an explicit roster order; no graph alias escapes. */
		static Loaded combine(final List<Loaded> families,final List<String> roster){
			if(families.size()!=roster.size() || families.isEmpty()){throw new IllegalArgumentException("Combined family count differs from the required roster");}
			final String[] ids=new String[roster.size()];final AAGraph[] combined=new AAGraph[roster.size()];
			final HashSet<String> seen=new HashSet<String>();
			for(int i=0; i<ids.length; i++){
				final Loaded family=families.get(i);ids[i]=roster.get(i);
				if(family.graphs.length!=1 || !ids[i].equals(family.repIds[0]) || !seen.add(ids[i])){
					throw new IllegalArgumentException("Missing, duplicate, or out-of-order combined family at "+i);
				}
				combined[i]=family.graphs[0];
			}
			return new Loaded(ids,combined);
		}
		/** Writes this immutable loaded model set without exposing it to the caller. */
		void writeBundle(final Path out,final byte[][] provenance) throws IOException{
			final java.util.ArrayList<HbmBundleBuilder.FamilyInput> families=new java.util.ArrayList<HbmBundleBuilder.FamilyInput>();
			for(int i=0; i<graphs.length; i++){families.add(new HbmBundleBuilder.FamilyInput(repIds[i],graphs[i].pivot,graphs[i]));}
			HbmBundleBuilder.build(out,families,provenance);
		}
		/** Uses an explicitly frozen background when comparing a rebuilt library with its parent. */
		HbmPositionModel[] positionModels(String kind,double beta,boolean clip,double clipMin,double[] background){
			final HbmPositionModel[] models=new HbmPositionModel[graphs.length];
			for(int i=0; i<graphs.length; i++){models[i]=HbmPositionModel.deriveOne(graphs[i],kind,beta,clip,background,clipMin);}
			return models;
		}
		/** Scores a query against family {@code i}'s private graph, applying the AAGraphScorer "final
		 *  conditional scorer" gate: returns {@code null} unless {@code passedAllOtherFilters} is true
		 *  (the caller owns that decision — this class evaluates no filter itself, matching
		 *  {@link AAGraphScorer#scoreIfSurvivor}'s own contract). Returns a fresh {@code {raw,relative}}
		 *  array when the gate passes — never a graph alias. */
		public float[] score(final int i, final byte[] queryEnc, final AAAlignment aln, final boolean passedAllOtherFilters){
			return AAGraphScorer.scoreIfSurvivor(graphs[i], queryEnc, aln, passedAllOtherFilters);
		}
		/** TEST ONLY (package-private): asserts family {@code i}'s private graph structurally equals {@code expected}
		 *  (knobs, every node field, and both REF and DEL ins chains) without exposing the graph. */
		void assertStructuralMatch(final int i, final AAGraph expected){ structuralEqual(graphs[i], expected); }
		/** Test-only exact comparison without exposing either bundle's private graphs. */
		void assertStructuralEquivalent(final Loaded expected){
			if(!Arrays.equals(repIds, expected.repIds)){throw new AssertionError("bundle family order differs");}
			for(int i=0; i<graphs.length; i++){structuralEqual(graphs[i], expected.graphs[i]);}
		}
	}

	/**
	 * Loads a bundle after checking either the eight runtime semantic hashes or
	 * all sixteen runtime-and-build provenance hashes.  Production assignment may
	 * supply the runtime subset; calibration and provenance-sensitive tools should
	 * supply all sixteen so the builder, roster, corpus, and membership lineage is
	 * also enforced.
	 */
	public static Loaded load(final Path file, final List<String> roster,
			final ConsensusProvider consensus, final byte[][] trustedProvenanceHashes)
			throws IOException{
		return load(file, roster, consensus, trustedProvenanceHashes, 1);
	}

	/**
	 * Loads a bundle and optionally suppresses REF/INS residue histogram slots
	 * whose counts are below {@code minCount}.  This is an in-memory experimental
	 * scorer perturbation; the on-disk bundle and its provenance remain hash-bound.
	 */
	public static Loaded load(final Path file, final List<String> roster,
			final ConsensusProvider consensus, final byte[][] trustedProvenanceHashes,
			final int minCount) throws IOException{
		if(minCount<1){throw new IllegalArgumentException("minCount must be at least 1: "+minCount);}
		if(HbmCompactTextReader.matches(file)){
			return HbmCompactTextReader.loadVerified(file.toString(), roster, consensus, trustedProvenanceHashes, minCount);
		}
		if(HbmDenseTextBundle.matches(file)){
			return HbmDenseTextLoader.loadVerified(file.toString(), roster, consensus, trustedProvenanceHashes,
				Math.max(1, Math.min(256, shared.Shared.threads())), minCount);
		}
		if(HbmA48Codec.isEncoded(file)){
			final Path temporary=Files.createTempDirectory("mqhb-a48-");
			final Path decoded=temporary.resolve("decoded.mqhb");
			try{
				HbmA48Codec.decode(file, decoded);
				return load(decoded, roster, consensus, trustedProvenanceHashes, minCount);
			}finally{Files.deleteIfExists(decoded); Files.deleteIfExists(temporary);}
		}
		//Roster validation BEFORE any allocation: nonempty, unique, strict-UTF-8 ids.
		if(roster==null || roster.isEmpty()){throw new IllegalArgumentException("nonempty roster required");}
		final HashSet<String> seen=new HashSet<String>();
		for(final String id : roster){
			if(id==null || id.isEmpty()){throw new IllegalArgumentException("blank roster id");}
			utf8StrictEncode(id);//rejects unpaired surrogates / unencodable
			if(!seen.add(id)){throw new IllegalArgumentException("duplicate roster id: "+id);}
		}
		if(consensus==null){throw new IllegalArgumentException("null consensus provider");}
		if(trustedProvenanceHashes==null ||
				(trustedProvenanceHashes.length!=HbmBundleFormat.RUNTIME_SEMANTIC_COUNT &&
				 trustedProvenanceHashes.length!=HbmBundleFormat.PROVENANCE_COUNT)){
			throw new IllegalArgumentException("trustedProvenanceHashes must have "+
				HbmBundleFormat.RUNTIME_SEMANTIC_COUNT+" or "+
				HbmBundleFormat.PROVENANCE_COUNT+" entries");
		}
		for(final byte[] tr : trustedProvenanceHashes){
			if(tr==null || (tr.length!=DigestSuffix.BYTE_LENGTH && tr.length!=HbmBundleFormat.SHA256_LEN)){
				throw new IllegalArgumentException("bad trusted hash length");
			}
		}

		final long fileLen=Files.size(file);
		if(fileLen<HbmBundleFormat.FIXED_HEADER_LEN){throw new IllegalArgumentException("file shorter than header");}

		try(FileChannel ch=FileChannel.open(file, StandardOpenOption.READ)){
			final byte[] header=readAt(ch, 0, HbmBundleFormat.FIXED_HEADER_LEN, fileLen);
			if(!HbmBundleFormat.magicMatches(header, HbmBundleFormat.OFF_MAGIC)){throw new IllegalArgumentException("bad magic");}
			if(u16(header, HbmBundleFormat.OFF_FORMAT_VERSION)!=HbmBundleFormat.FORMAT_VERSION){throw new IllegalArgumentException("unsupported format_version");}
			if(u16(header, HbmBundleFormat.OFF_NAA)!=HbmBundleFormat.NAA){throw new IllegalArgumentException("naa mismatch");}
			if(u16(header, HbmBundleFormat.OFF_PAD)!=HbmBundleFormat.EXPECT_PAD){throw new IllegalArgumentException("pad!=0");}
			final int n=i32(header, HbmBundleFormat.OFF_N_FAMILIES);
			if(n<1){throw new IllegalArgumentException("n_families<1");}
			if(n!=roster.size()){throw new IllegalArgumentException("n_families "+n+" != roster "+roster.size());}
			//Bound n_families BEFORE allocating arrays: reserve minimum directory bytes AND one minimum
			//L=1 family block per family AND the 36-byte trailer (v4 — v3 only reserved directory bytes).
			final long minPerFamily=(HbmBundleFormat.DIR_ENTRY_FIXED_BYTES+1)+MIN_BLOCK_BYTES;
			final long trailerBytes=HbmBundleFormat.SHA256_LEN+HbmBundleFormat.MAGIC_LEN;
			if((long)n*minPerFamily+trailerBytes>fileLen-HbmBundleFormat.FIXED_HEADER_LEN){throw new IllegalArgumentException("n_families too large for file size");}
			checkKnobs(header);
			//§4: runtime-semantic compatibility FIRST, before any graph.
			for(int k=0; k<trustedProvenanceHashes.length; k++){
				final int off=HbmBundleFormat.OFF_PROVENANCE+k*HbmBundleFormat.SHA256_LEN;
				final byte[] trusted=trustedProvenanceHashes[k];
				final int trustedOffset=HbmBundleFormat.SHA256_LEN-trusted.length;
				for(int b=0; b<trusted.length; b++){
					if(header[off+trustedOffset+b]!=trusted[b]){
						throw new IllegalArgumentException(
							"trusted provenance hash mismatch at provenance["+k+"]");
					}
				}
			}
			if(HbmBundleFormat.isAllZero(header, HbmBundleFormat.OFF_PROVENANCE+HbmBundleFormat.PROV_LOADER_SOURCE*HbmBundleFormat.SHA256_LEN)
					|| HbmBundleFormat.isAllZero(header, HbmBundleFormat.OFF_PROVENANCE+HbmBundleFormat.PROV_LOADER_CLASS*HbmBundleFormat.SHA256_LEN)){
				throw new IllegalArgumentException("all-zero loader_source/class hash in a v1 bundle");
			}

			//Directory.
			final long[] blockOffset=new long[n], blockLen=new long[n], blockCrc=new long[n];
			long pos=HbmBundleFormat.DIRECTORY_START;
			for(int i=0; i<n; i++){
				final int repLen=i32(readAt(ch, pos, 4, fileLen), 0);
				if(repLen<1 || repLen>HbmBundleFormat.CAP_REPID){throw new IllegalArgumentException("rep_id_len out of range at "+i);}
				pos+=4;
				final byte[] repBuf=readAt(ch, pos, repLen, fileLen); pos+=repLen;
				final String repId=utf8StrictDecode(repBuf);//strict: malformed rep bytes → reject
				if(!repId.equals(roster.get(i))){throw new IllegalArgumentException("directory rep_id != roster at rank "+i);}
				final byte[] tri=readAt(ch, pos, HbmBundleFormat.DIR_ENTRY_FIXED_BYTES-4, fileLen); pos+=HbmBundleFormat.DIR_ENTRY_FIXED_BYTES-4;
				blockOffset[i]=i64(tri, 0); blockLen[i]=i64(tri, 8); blockCrc[i]=u32(tri, 16);
			}
			final long directoryEnd=pos;
			//Gap-free layout + overflow-safe bounds + block-len cap + CAP_TOTAL_BYTES (Σblock_len only).
			long expect=directoryEnd, sumBlockBytes=0;
			for(int i=0; i<n; i++){
				if(blockOffset[i]!=expect){throw new IllegalArgumentException("block "+i+" not gap-free");}
				if(blockLen[i]<1 || blockLen[i]>CAP_BLOCK_BYTES){throw new IllegalArgumentException("block "+i+" length out of [1,CAP_BLOCK_BYTES]");}
				if(blockOffset[i]<0 || blockOffset[i]>fileLen || blockLen[i]>fileLen-blockOffset[i]){throw new IllegalArgumentException("block "+i+" out of bounds");}
				if(blockLen[i]>HbmBundleFormat.CAP_TOTAL_BYTES-sumBlockBytes){throw new IllegalArgumentException("Σblock_len exceeds CAP_TOTAL_BYTES");}
				sumBlockBytes+=blockLen[i];
				expect+=blockLen[i];
			}
			final long trailerStart=expect;
			if(fileLen!=trailerStart+HbmBundleFormat.SHA256_LEN+HbmBundleFormat.MAGIC_LEN){throw new IllegalArgumentException("trailing bytes or short file");}
			//Trailer: prefix_sha256 + sentinel.
			final byte[] trailer=readAt(ch, trailerStart, HbmBundleFormat.SHA256_LEN+HbmBundleFormat.MAGIC_LEN, fileLen);
			final byte[] prefixSha=streamingSha(ch, 0, trailerStart);
			for(int b=0; b<HbmBundleFormat.SHA256_LEN; b++){if(trailer[b]!=prefixSha[b]){throw new IllegalArgumentException("prefix_sha256 mismatch");}}
			if(!HbmBundleFormat.magicMatches(trailer, HbmBundleFormat.SHA256_LEN)){throw new IllegalArgumentException("bad end sentinel");}

			//Blocks — eager reconstruction into private graphs.
			final AAGraph[] graphs=new AAGraph[n];
			final String[] repIds=new String[n];
			long totalNodes=0;
			for(int i=0; i<n; i++){
				final byte[] block=readAt(ch, blockOffset[i], (int)blockLen[i], fileLen);//blockLen<=CAP_BLOCK_BYTES<INT_MAX
				if(HbmBundleFormat.crc32(block, 0, block.length)!=blockCrc[i]){throw new IllegalArgumentException("block "+i+" CRC32 mismatch");}
				totalNodes=reconstruct(block, roster.get(i), consensus, graphs, i, totalNodes, minCount);
				repIds[i]=roster.get(i);
			}
			return new Loaded(repIds, graphs);
		}
	}

	private static long reconstruct(final byte[] block, final String repId, final ConsensusProvider consensus,
			final AAGraph[] graphs, final int gi, long totalNodes, final int minCount){
		final Cursor c=new Cursor(block);
		final int L=c.i32();
		if(L<1 || L>HbmBundleFormat.CAP_L){throw new IllegalArgumentException("L out of range: "+L);}
		final byte[] cSha=c.bytes(HbmBundleFormat.SHA256_LEN);
		//Clone the provider's consensus FIRST, then hash the clone, then build on the clone (no provider alias).
		final byte[] provided=consensus.consensusFor(repId);
		if(provided==null || provided.length!=L){throw new IllegalArgumentException("consensus for '"+repId+"' missing/wrong length");}
		final byte[] clone=Arrays.copyOf(provided, provided.length);
		if(!Arrays.equals(cSha, HbmBundleFormat.sha256(clone, 0, clone.length))){throw new IllegalArgumentException("consensus_sha256 mismatch for '"+repId+"'");}
		//Every residue must be valid BEFORE AAGraph construction (v4): stop (20) would otherwise build
		//silently, and negative/above-X_CODE bytes would otherwise crash AAGraph's own scaffold with an
		//unexplained ArrayIndexOutOfBoundsException instead of this explicit, diagnosable rejection.
		for(int k=0; k<L; k++){if(!isValidResidue(clone[k])){throw new IllegalArgumentException("invalid consensus residue for '"+repId+"' at "+k+": "+clone[k]);}}
		//Cap the 2*L REF+DEL nodes BEFORE allocating the graph.
		if(2L*L>HbmBundleFormat.CAP_TOTAL_NODES-totalNodes){throw new IllegalArgumentException("CAP_TOTAL_NODES exceeded");}
		final AAGraph g=new AAGraph(clone, 0);
		checkGraphKnobs(g);
		totalNodes+=2L*L;
		for(int i=0; i<L; i++){
			final int refCountSum=c.i32();
			if(refCountSum<1){throw new IllegalArgumentException("REF countSum<1 at "+i);}
			overwrite(g.ref[i], c, refCountSum, minCount, clone[i]&0xff);
			if(g.ref[i].count[clone[i]&0xff]<1){throw new IllegalArgumentException("REF pivot count<1 at "+i);}
			totalNodes=readChain(c, g.ref[i], i+1, totalNodes, minCount);
			final int delCountSum=c.i32();
			if(delCountSum<0){throw new IllegalArgumentException("DEL countSum<0 at "+i);}
			g.del[i].countSum=delCountSum; g.del[i].weightSum=delCountSum;
			totalNodes=readChain(c, g.del[i], i+1, totalNodes, minCount);
		}
		c.end();//block must be fully consumed (no trailing block bytes)
		graphs[gi]=g;
		return totalNodes;
	}

	/** True iff {@code b} is a valid encoded protein residue for a consensus byte: real residues 0-19,
	 *  or the unknown code {@link Blosum62#X_CODE} (21). Duplicated identically in
	 *  {@code HbmBundleBuilder} — see that class's javadoc for why it is not promoted into the
	 *  already-accepted-and-applied {@code HbmBundleFormat} (Increment 1). */
	private static boolean isValidResidue(final byte b){
		return (b>=0 && b<=19) || b==Blosum62.X_CODE;
	}

	/** Overwrites a REF/INS node's count/weight/countSum/weightSum from the stored values; Σcount==countSum. */
	private static void overwrite(final AAGraphNode node, final Cursor c,
			final int countSum, final int minCount, final int protectedResidue){
		long sum=0;
		int retainedSum=0;
		for(int k=0; k<HbmBundleFormat.NAA; k++){
			final int stored=c.i32();
			if(stored<0){throw new IllegalArgumentException("negative count");}
			final int retained=(stored>=minCount || k==protectedResidue) ? stored : 0;
			node.count[k]=retained; node.weight[k]=retained;
			sum+=stored;
			retainedSum=Math.addExact(retainedSum, retained);
		}
		if(sum!=countSum){throw new IllegalArgumentException("countSum "+countSum+" != Σcount "+sum);}
		node.countSum=retainedSum; node.weightSum=retainedSum;
	}

	/** Reads an ins chain off {@code parent}; checks the total-node cap AND the declared chain length
	 *  against the complete remaining block bytes BEFORE allocating each INS node (v4: the remaining-
	 *  bytes check used to be implicit in the per-read {@code Cursor.need} bounds check, which could
	 *  fire only AFTER earlier nodes in the same chain were already allocated). */
	private static long readChain(final Cursor c, final AAGraphNode parent,
			final int rpos, long totalNodes, final int minCount){
		final int chainLen=c.i32();
		if(chainLen<0 || chainLen>HbmBundleFormat.CAP_CHAIN){throw new IllegalArgumentException("chain_len out of range: "+chainLen);}
		final long neededBytes=(long)chainLen*HbmBundleFormat.INS_RECORD_BYTES;
		if(neededBytes>c.remaining()){throw new IllegalArgumentException("chain_len "+chainLen+" exceeds remaining block bytes");}
		AAGraphNode prev=null;
		boolean truncated=false;
		for(int j=0; j<chainLen; j++){
			final int countSum=c.i32();
			if(countSum<1){throw new IllegalArgumentException("INS countSum<1");}
			if(1L>HbmBundleFormat.CAP_TOTAL_NODES-totalNodes){throw new IllegalArgumentException("CAP_TOTAL_NODES exceeded");}
			final AAGraphNode node=new AAGraphNode(Blosum62.X_CODE, AAGraphNode.INS, rpos);
			totalNodes++;
			overwrite(node, c, countSum, minCount, -1);
			if(node.countSum<1){truncated=true;}
			else if(!truncated){
				if(prev==null){parent.insEdge=node;}else{prev.insEdge=node;}
				prev=node;
			}
		}
		return totalNodes;
	}

	/*--------------------------------------------------------------*/
	/*----------------   Structural compare (test)   ----------------*/
	/*--------------------------------------------------------------*/

	private static void structuralEqual(final AAGraph a, final AAGraph b){
		if(a.pad!=b.pad || a.ref.length!=b.ref.length || !Arrays.equals(a.pivot, b.pivot)){throw new AssertionError("graph pad/L/pivot differ");}
		if(a.weightByIdentity!=b.weightByIdentity || a.baseWeight!=b.baseWeight || a.minDepth!=b.minDepth
				|| fb(a.MAF_sub)!=fb(b.MAF_sub) || fb(a.MAF_del)!=fb(b.MAF_del) || fb(a.MAF_ins)!=fb(b.MAF_ins)
				|| fb(a.identityCeiling)!=fb(b.identityCeiling) || fb(a.trimDepthFraction)!=fb(b.trimDepthFraction)){throw new AssertionError("graph knobs differ");}
		for(int i=0; i<a.ref.length; i++){
			nodeEqual(a.ref[i], b.ref[i]);
			chainEqual(a.ref[i].insEdge, b.ref[i].insEdge);
			nodeEqual(a.del[i], b.del[i]);
			chainEqual(a.del[i].insEdge, b.del[i].insEdge);
		}
	}
	private static void chainEqual(AAGraphNode x, AAGraphNode y){
		while(x!=null && y!=null){nodeEqual(x, y); x=x.insEdge; y=y.insEdge;}
		if(x!=null || y!=null){throw new AssertionError("ins chain length differs");}
	}
	private static void nodeEqual(final AAGraphNode a, final AAGraphNode b){
		if(a.type!=b.type || a.rpos!=b.rpos || a.refResidue!=b.refResidue ||
				a.countSum!=b.countSum || a.weightSum!=b.weightSum){
			throw new AssertionError("node meta differs: actual type="+a.type+
				" rpos="+a.rpos+" ref="+a.refResidue+" countSum="+a.countSum+
				" weightSum="+a.weightSum+" expected type="+b.type+" rpos="+
				b.rpos+" ref="+b.refResidue+" countSum="+b.countSum+
				" weightSum="+b.weightSum);
		}
		if((a.count==null)!=(b.count==null) || (a.weight==null)!=(b.weight==null)){throw new AssertionError("node array nullity differs");}
		if(a.count!=null){for(int k=0; k<HbmBundleFormat.NAA; k++){if(a.count[k]!=b.count[k] || a.weight[k]!=b.weight[k]){throw new AssertionError("node counts differ");}}}
	}
	private static int fb(final float f){return Float.floatToRawIntBits(f);}

	static void checkGraphKnobs(final AAGraph g){
		if(g.pad!=HbmBundleFormat.EXPECT_PAD || g.weightByIdentity!=HbmBundleFormat.KNOB_WEIGHT_BY_IDENTITY
				|| g.baseWeight!=HbmBundleFormat.KNOB_BASE_WEIGHT || g.minDepth!=HbmBundleFormat.KNOB_MIN_DEPTH
				|| fb(g.MAF_sub)!=fb(HbmBundleFormat.KNOB_MAF_SUB) || fb(g.MAF_del)!=fb(HbmBundleFormat.KNOB_MAF_DEL)
				|| fb(g.MAF_ins)!=fb(HbmBundleFormat.KNOB_MAF_INS) || fb(g.identityCeiling)!=fb(HbmBundleFormat.KNOB_IDENTITY_CEILING)
				|| fb(g.trimDepthFraction)!=fb(HbmBundleFormat.KNOB_TRIM_DEPTH_FRACTION)){
			throw new IllegalArgumentException("reconstructed graph knobs are not canonical");
		}
	}
	private static void checkKnobs(final byte[] h){
		if(h[HbmBundleFormat.OFF_WEIGHT_BY_IDENTITY]!=(HbmBundleFormat.KNOB_WEIGHT_BY_IDENTITY?1:0)
				|| i32(h, HbmBundleFormat.OFF_BASE_WEIGHT)!=HbmBundleFormat.KNOB_BASE_WEIGHT
				|| i32(h, HbmBundleFormat.OFF_MAF_SUB)!=fb(HbmBundleFormat.KNOB_MAF_SUB)
				|| i32(h, HbmBundleFormat.OFF_MAF_DEL)!=fb(HbmBundleFormat.KNOB_MAF_DEL)
				|| i32(h, HbmBundleFormat.OFF_MAF_INS)!=fb(HbmBundleFormat.KNOB_MAF_INS)
				|| i32(h, HbmBundleFormat.OFF_MIN_DEPTH)!=HbmBundleFormat.KNOB_MIN_DEPTH
				|| i32(h, HbmBundleFormat.OFF_IDENTITY_CEILING)!=fb(HbmBundleFormat.KNOB_IDENTITY_CEILING)
				|| i32(h, HbmBundleFormat.OFF_TRIM_DEPTH_FRACTION)!=fb(HbmBundleFormat.KNOB_TRIM_DEPTH_FRACTION)){
			throw new IllegalArgumentException("header knobs are not canonical");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------   Bounded byte-cursor + I/O   ----------------*/
	/*--------------------------------------------------------------*/

	/** Cursor over a block byte[] with overflow-safe remaining-byte checks on every read. */
	static final class Cursor{
		final byte[] a; int p;
		Cursor(final byte[] a){this.a=a;}
		private void need(final int len){if(len<0 || len>a.length-p){throw new IllegalArgumentException("read past end of block (need "+len+" at "+p+"/"+a.length+")");}}
		int i32(){need(4); final int v=((a[p]&0xff)<<24)|((a[p+1]&0xff)<<16)|((a[p+2]&0xff)<<8)|(a[p+3]&0xff); p+=4; return v;}
		byte[] bytes(final int len){need(len); final byte[] out=Arrays.copyOfRange(a, p, p+len); p+=len; return out;}
		/** Bytes remaining after the cursor's current position (always &ge; 0). */
		int remaining(){return a.length-p;}
		void end(){if(p!=a.length){throw new IllegalArgumentException("block has trailing bytes (consumed "+p+" of "+a.length+")");}}
	}

	private static byte[] readAt(final FileChannel ch, final long pos, final int len, final long fileLen) throws IOException{
		if(pos<0 || len<0 || pos>fileLen || len>fileLen-pos){throw new IllegalArgumentException("read out of bounds pos="+pos+" len="+len);}
		final ByteBuffer bb=ByteBuffer.allocate(len); long at=pos;
		while(bb.hasRemaining()){final int r=ch.read(bb, at); if(r<0){throw new IOException("EOF at "+at);} at+=r;}
		return bb.array();
	}
	private static byte[] streamingSha(final FileChannel ch, final long from, final long to) throws IOException{
		final MessageDigest md; try{md=MessageDigest.getInstance("SHA-256");}catch(NoSuchAlgorithmException e){throw new RuntimeException("SHA-256 unavailable", e);}
		final ByteBuffer buf=ByteBuffer.allocate(1<<16); long pos=from;
		while(pos<to){buf.clear(); final long rem=to-pos; if(rem<buf.capacity()){buf.limit((int)rem);}
			final int r=ch.read(buf, pos); if(r<0){throw new IOException("EOF hashing prefix at "+pos);} buf.flip(); md.update(buf); pos+=r;}
		return md.digest();
	}
	private static int u16(final byte[] a, final int o){return ((a[o]&0xff)<<8)|(a[o+1]&0xff);}
	private static int i32(final byte[] a, final int o){return ((a[o]&0xff)<<24)|((a[o+1]&0xff)<<16)|((a[o+2]&0xff)<<8)|(a[o+3]&0xff);}
	private static long u32(final byte[] a, final int o){return i32(a, o)&0xffffffffL;}
	private static long i64(final byte[] a, final int o){long v=0; for(int k=0; k<8; k++){v=(v<<8)|(a[o+k]&0xff);} return v;}

	//Strict UTF-8 (same approach as the accepted HbmMemberManifest).
	static byte[] utf8StrictEncode(final String s){
		final CharsetEncoder enc=StandardCharsets.UTF_8.newEncoder().onMalformedInput(CodingErrorAction.REPORT).onUnmappableCharacter(CodingErrorAction.REPORT);
		try{final ByteBuffer bb=enc.encode(CharBuffer.wrap(s)); final byte[] out=new byte[bb.remaining()]; bb.get(out); return out;}
		catch(CharacterCodingException e){throw new IllegalArgumentException("id not valid UTF-16/UTF-8: "+e.getMessage());}
	}
	static String utf8StrictDecode(final byte[] raw){
		final CharsetDecoder dec=StandardCharsets.UTF_8.newDecoder().onMalformedInput(CodingErrorAction.REPORT).onUnmappableCharacter(CodingErrorAction.REPORT);
		try{return dec.decode(ByteBuffer.wrap(raw)).toString();}
		catch(CharacterCodingException e){throw new IllegalArgumentException("malformed UTF-8 rep_id: "+e.getMessage());}
	}
}
