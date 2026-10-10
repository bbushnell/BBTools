package prot;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import parse.Parser;
import structures.ByteBuilder;

/**
 * Produces a separate MQHBv1 rare-state candidate under the D242 contract.
 * Rarity is strictly below the selected cutoff of original scorer denominators;
 * cutoff=0.005 is the default and cutoff=0.01 is the only alternate. Pivot
 * support and the native insertion fallback remain intact. No production model
 * is edited.
 * @author Nilou
 */
public final class HbmRareStateProbe {

	/** Parses named options; fresh output directories preserve the accepted input. */
	public static void main(final String[] args) throws Exception{
		String input=null, reference=null, provenance=null, output=null, mode="prune", cutoff="0.005";
		boolean modeSeen=false, cutoffSeen=false;
		final String[] expanded=Parser.parseConfig(args);
		if(expanded.length==1 && expanded[0].equalsIgnoreCase("selftest=t")){
			HbmRareStateProbeTest.main(new String[0]); return;
		}
		for(String arg:expanded){
			final int equals=arg.indexOf('=');
			if(equals<1 || equals==arg.length()-1){throw new IllegalArgumentException("Expected flag=value: "+arg);}
			final String key=arg.substring(0, equals), value=arg.substring(equals+1);
			if(key.equalsIgnoreCase("in") && input==null){input=value;}
			else if(key.equalsIgnoreCase("ref") && reference==null){reference=value;}
			else if(key.equalsIgnoreCase("provenance") && provenance==null){provenance=value;}
			else if(key.equalsIgnoreCase("out") && output==null){output=value;}
			else if(key.equalsIgnoreCase("mode") && !modeSeen){mode=value; modeSeen=true;}
			else if(key.equalsIgnoreCase("cutoff") && !cutoffSeen){cutoff=value; cutoffSeen=true;}
			else{throw new IllegalArgumentException("Unknown or duplicate option: "+key);}
		}
		if(input==null || reference==null || provenance==null || output==null ||
				!(mode.equals("prune") || mode.equals("identity"))){
			throw new IllegalArgumentException("Required in= ref= provenance= out=fresh-directory [mode=prune|identity] [cutoff=0.005|0.01]");
		}
		final int cutoffScale=parseCutoff(cutoff);
		run(Paths.get(input), reference, provenance, Paths.get(output), mode.equals("prune"), cutoffScale);
	}

	private HbmRareStateProbe(){}

	/** Native-loads both artifacts; identity mode verifies deterministic reconstruction. */
	static void run(final Path input, final String reference, final String provenancePath,
			final Path output, final boolean remove) throws Exception{
		run(input, reference, provenancePath, output, remove, 200);
	}

	/** Runs with an explicit supported cutoff scale: 200 for0.5%,100 for1%. */
	static void run(final Path input, final String reference, final String provenancePath,
			final Path output, final boolean remove, final int cutoffScale) throws Exception{
		validateCutoff(cutoffScale);
		if(Files.exists(output)){throw new IllegalArgumentException("Fresh output directory required: "+output);}
		final byte[] originalDigest=digest(input);
		final List<ProteinSequence> sequences=ProteinSearch.readFasta(reference);
		final ArrayList<String> roster=new ArrayList<String>(sequences.size());
		final HashMap<String, byte[]> consensus=new HashMap<String, byte[]>(sequences.size()*2);
		for(ProteinSequence sequence:sequences){
			roster.add(sequence.id);
			if(consensus.put(sequence.id, sequence.enc)!=null){throw new IllegalArgumentException("Duplicate consensus: "+sequence.id);}
		}
		final byte[][] trusted=HbmBundleLoader.loadSemanticProvenance(provenancePath);
		validateLoad(input, roster, consensus, trusted, null);
		final Stats totals=new Stats();
		final ByteBuilder report=new ByteBuilder(1<<20);
		report.append("family\tresidue_slots_removed\tresidue_count_removed\tpivot_slots_protected\tdel_states_removed\tdel_count_removed\tins_nodes_removed\tins_count_removed\n");
		final ArrayList<HbmBundleBuilder.FamilyInput> families=new ArrayList<HbmBundleBuilder.FamilyInput>(roster.size());
		final byte[][] sourceProvenance=new byte[HbmBundleFormat.PROVENANCE_COUNT][];
		try(FileChannel channel=FileChannel.open(input, StandardOpenOption.READ)){
			final byte[] header=read(channel, 0, HbmBundleFormat.FIXED_HEADER_LEN);
			for(int i=0; i<sourceProvenance.length; i++){
				final int start=HbmBundleFormat.OFF_PROVENANCE+i*HbmBundleFormat.SHA256_LEN;
				sourceProvenance[i]=Arrays.copyOfRange(header, start, start+HbmBundleFormat.SHA256_LEN);
			}
			long position=HbmBundleFormat.DIRECTORY_START;
			for(int rank=0; rank<roster.size(); rank++){
				final int nameLength=ByteBuffer.wrap(read(channel, position, 4)).getInt(); position+=4;
				if(nameLength<1 || nameLength>HbmBundleFormat.CAP_REPID){throw new IOException("Invalid representative length");}
				final String name=new String(read(channel, position, nameLength), StandardCharsets.UTF_8); position+=nameLength;
				if(!name.equals(roster.get(rank))){throw new IOException("Consensus/directory order differs at "+rank);}
				final ByteBuffer entry=ByteBuffer.wrap(read(channel, position, 20)); position+=20;
				final long offset=entry.getLong(), length=entry.getLong(), crc=entry.getInt()&0xffffffffL;
				if(length<36 || length>HbmBundleLoader.CAP_BLOCK_BYTES){throw new IOException("Invalid block size for "+name);}
				final byte[] block=read(channel, offset, (int)length);
				if(HbmBundleFormat.crc32(block, 0, block.length)!=crc){throw new IOException("Changed block CRC for "+name);}
				final AAGraph graph=reconstruct(block, consensus.get(name));
				final HbmBundleBuilder.FamilyInput family=new HbmBundleBuilder.FamilyInput(name, consensus.get(name), graph);
				HbmBundleBuilder.validateGraph(family);
				final Stats stats=remove ? prune(graph, cutoffScale) : new Stats();
				HbmBundleBuilder.validateGraph(family);
				families.add(family); totals.add(stats); stats.append(report, name);
			}
		}
		if(!Arrays.equals(originalDigest, digest(input))){throw new IOException("Input changed while reconstructing");}
		Files.createDirectory(output);
		final Path candidate=output.resolve("candidate.mqhb");
		HbmBundleBuilder.build(candidate, families, sourceProvenance);
		validateLoad(candidate, roster, consensus, trusted, families);
		if(!remove && !Arrays.equals(originalDigest, digest(candidate))){
			throw new IOException("Identity reconstruction differs; do not interpret pruning results");
		}
		if(!Arrays.equals(originalDigest, digest(input))){throw new IOException("Original input changed");}
		Files.write(output.resolve("families.tsv"), report.toBytes(), StandardOpenOption.CREATE_NEW);
		final ByteBuilder summary=new ByteBuilder(4096);
		summary.append("metric\tvalue\nmode\t").append(remove ? "prune" : "identity").nl();
		summary.append("cutoff\t").append(cutoffScale==200 ? "0.005" : "0.01").nl();
		summary.append("cutoff_denominator\t").append(cutoffScale).nl();
		summary.append("families\t").append(families.size()).nl();
		summary.append("input_bytes\t").append(Files.size(input)).nl();
		summary.append("output_bytes\t").append(Files.size(candidate)).nl();
		summary.append("input_sha80\t").append(sha80(originalDigest)).nl();
		summary.append("output_sha80\t").append(sha80(digest(candidate))).nl();
		summary.append("residue_slots_removed\t").append(totals.residueSlots).nl();
		summary.append("residue_count_removed\t").append(totals.residueMass).nl();
		summary.append("pivot_slots_protected\t").append(totals.protectedPivots).nl();
		summary.append("del_states_removed\t").append(totals.delStates).nl();
		summary.append("del_count_removed\t").append(totals.delMass).nl();
		summary.append("ins_nodes_removed\t").append(totals.insNodes).nl();
		summary.append("ins_count_removed\t").append(totals.insMass).nl();
		summary.append("status\tNATIVE_LOAD_AND_GRAPH_COMPARE_PASS\n");
		Files.write(output.resolve("summary.tsv"), summary.toBytes(), StandardOpenOption.CREATE_NEW);
		System.out.print(summary.toString());
	}

	/** Uses the production loader, including its private graph structural comparator. */
	private static void validateLoad(final Path path, final List<String> roster,
			final HashMap<String, byte[]> consensus, final byte[][] trusted,
			final List<HbmBundleBuilder.FamilyInput> expected) throws IOException{
		final HbmBundleLoader.Loaded loaded=HbmBundleLoader.load(path, roster, consensus::get, trusted);
		if(loaded.familyCount()!=roster.size()){throw new IllegalStateException("Native loader family count changed");}
		if(expected!=null){for(int i=0; i<expected.size(); i++){loaded.assertStructuralMatch(i, expected.get(i).graph);}}
	}

	/** Reconstructs an already native-validated block without exposing loader-owned graphs. */
	static AAGraph reconstruct(final byte[] bytes, final byte[] consensus) throws IOException{
		final ByteBuffer block=ByteBuffer.wrap(bytes);
		final int length=block.getInt();
		if(length!=consensus.length){throw new IOException("Consensus length changed");}
		final byte[] expected=new byte[32]; block.get(expected);
		if(!Arrays.equals(expected, HbmBundleFormat.sha256(consensus, 0, length))){throw new IOException("Consensus digest differs");}
		final AAGraph graph=new AAGraph(consensus.clone(), 0);
		for(int i=0; i<length; i++){
			readHistogram(block, graph.ref[i]); readChain(block, graph.ref[i], i+1);
			graph.del[i].countSum=graph.del[i].weightSum=block.getInt();
			readChain(block, graph.del[i], i+1);
		}
		if(block.hasRemaining()){throw new IOException("Unconsumed family block bytes");}
		return graph;
	}

	/** Restores both integer arrays; the native builder checks sums and positivity. */
	private static void readHistogram(final ByteBuffer block, final AAGraphNode node){
		node.countSum=node.weightSum=block.getInt();
		for(int k=0; k<node.count.length; k++){node.count[k]=node.weight[k]=block.getInt();}
	}

	/** Checks the byte bound before allocating any insertion nodes. */
	private static void readChain(final ByteBuffer block, final AAGraphNode parent, final int position) throws IOException{
		final int length=block.getInt();
		if(length<0 || length>HbmBundleFormat.CAP_CHAIN || (long)length*HbmBundleFormat.INS_RECORD_BYTES>block.remaining()){
			throw new IOException("Invalid insertion chain length");
		}
		AAGraphNode tail=parent;
		for(int i=0; i<length; i++){
			final AAGraphNode node=new AAGraphNode(Blosum62.X_CODE, AAGraphNode.INS, position);
			readHistogram(block, node); tail.insEdge=node; tail=node;
		}
	}

	/** Applies the accepted strict cutoff using pre-mutation positional depths. */
	static Stats prune(final AAGraph graph){
		return prune(graph, 200);
	}

	/** Applies the selected strict cutoff to immutable pre-mutation denominators. */
	static Stats prune(final AAGraph graph, final int cutoffScale){
		validateCutoff(cutoffScale);
		final long[] depth=new long[graph.ref.length];
		for(int i=0; i<depth.length; i++){depth[i]=(long)graph.ref[i].countSum+graph.del[i].countSum;}
		final Stats stats=new Stats();
		for(int i=0; i<depth.length; i++){
			pruneChain(graph.ref[i], depth, stats, cutoffScale); pruneChain(graph.del[i], depth, stats, cutoffScale);
			pruneHistogram(graph.ref[i], graph.pivot[i], stats, cutoffScale);
			final AAGraphNode del=graph.del[i];
			if(rare(del.countSum, depth[i], cutoffScale)){
				stats.delStates++; stats.delMass+=del.countSum; del.countSum=del.weightSum=0;
			}
		}
		return stats;
	}

	/** Removes only a rare suffix; deeper insertions cannot survive by shifting left. */
	private static void pruneChain(final AAGraphNode parent, final long[] depth, final Stats stats, final int cutoffScale){
		int previousCount=Integer.MAX_VALUE;
		for(AAGraphNode node=parent.insEdge; node!=null; node=node.insEdge){
			if(node.rpos!=parent.rpos+1 || node.rpos>=depth.length || depth[node.rpos]<=0 ||
					node.countSum<1 || node.countSum>previousCount){
				throw new IllegalArgumentException("Insertion chain violates AAGraph.add prefix counts or scorer rpos denominator");
			}
			previousCount=node.countSum;
		}
		AAGraphNode previous=parent;
		for(AAGraphNode node=parent.insEdge; node!=null; node=node.insEdge){
			if(rare(node.countSum, depth[node.rpos], cutoffScale)){
				for(AAGraphNode removed=node; removed!=null; removed=removed.insEdge){
					stats.insNodes++; stats.insMass+=removed.countSum;
				}
				previous.insEdge=null; return;
			}
			pruneHistogram(node, -1, stats, cutoffScale); previous=node;
		}
	}

	/** Real-residue probabilities exclude X; protected REF pivots keep their full count. */
	private static void pruneHistogram(final AAGraphNode node, final int protectedResidue, final Stats stats, final int cutoffScale){
		long real=0;
		for(int k=0; k<20; k++){real+=node.count[k];}
		for(int k=0; k<20; k++){
			if(!rare(node.count[k], real, cutoffScale)){continue;}
			if(k==protectedResidue){stats.protectedPivots++; continue;}
			stats.residueSlots++; stats.residueMass+=node.count[k]; node.count[k]=node.weight[k]=0;
		}
		int sum=0;
		for(int count:node.count){sum=Math.addExact(sum, count);}
		node.countSum=node.weightSum=sum;
		assert(sum>0) : "selected cutoff removes only rare residues; native REF/INS nodes must retain support";
	}

	/** Default method compatibility: exact boundary and wide arithmetic at0.5%. */
	static boolean rare(final int count, final long denominator){return rare(count, denominator, 200);}

	/** Exact boundary and wide arithmetic for the supported cutoff scale. */
	static boolean rare(final int count, final long denominator, final int cutoffScale){
		validateCutoff(cutoffScale);
		return count>0 && (long)cutoffScale*count<denominator;
	}

	/** Parses only the two reviewed decimal cutoff spellings. */
	private static int parseCutoff(final String cutoff){
		if(cutoff.equals("0.005")){return 200;}
		if(cutoff.equals("0.01")){return 100;}
		throw new IllegalArgumentException("cutoff must be exactly0.005 or0.01: "+cutoff);
	}

	/** Prevents an accidental third cutoff from changing experiment semantics. */
	private static void validateCutoff(final int cutoffScale){
		if(cutoffScale!=200 && cutoffScale!=100){throw new IllegalArgumentException("cutoff scale must be200 or100");}
	}

	/** Positional binary I/O keeps directory offsets long and checks truncation. */
	private static byte[] read(final FileChannel channel, final long position, final int length) throws IOException{
		if(position<0 || length<0 || position>channel.size()-length){throw new IOException("Bundle read out of bounds");}
		final ByteBuffer buffer=ByteBuffer.allocate(length);
		while(buffer.hasRemaining()){
			if(channel.read(buffer, position+buffer.position())<=0){throw new IOException("No progress reading bundle");}
		}
		return buffer.array();
	}

	/** Hashes streamed bytes; callers never log the full SHA256 digest. */
	static byte[] digest(final Path path) throws Exception{
		final MessageDigest digest=MessageDigest.getInstance("SHA-256");
		try(FileChannel channel=FileChannel.open(path, StandardOpenOption.READ)){
			final ByteBuffer buffer=ByteBuffer.allocate(1<<20);
			while(channel.read(buffer)!=-1){buffer.flip(); digest.update(buffer); buffer.clear();}
		}
		return digest.digest();
	}

	/** Only the final80 bits are rendered, per the project's artifact-pin convention. */
	private static String sha80(final byte[] digest){
		final String hex=HbmMemberIndexFormat.toHexLower(digest);
		return hex.substring(hex.length()-20);
	}

	/** Separate event/mass counts avoid confusing removed observations with nodes. */
	static final class Stats {
		long residueSlots, residueMass, protectedPivots, delStates, delMass, insNodes, insMass;
		void add(final Stats other){
			residueSlots+=other.residueSlots; residueMass+=other.residueMass; protectedPivots+=other.protectedPivots;
			delStates+=other.delStates; delMass+=other.delMass; insNodes+=other.insNodes; insMass+=other.insMass;
		}
		void append(final ByteBuilder out, final String family){
			out.append(family).tab().append(residueSlots).tab().append(residueMass).tab().append(protectedPivots).tab()
				.append(delStates).tab().append(delMass).tab().append(insNodes).tab().append(insMass).nl();
		}
	}
}
