package prot.cluster;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashSet;
import java.util.Locale;
import java.util.Random;

import fileIO.ByteStreamWriter;
import prot.ProteinSequence;
import structures.ByteBuilder;

/** Writes near-equal contiguous shards of a seeded all-input permutation.
 * Parent ordinals, not text IDs, distinguish records with duplicate identifiers.
 * The caller must load with preserveFasta=true and reject aliases to source files;
 * this API has the loaded records, not the original source pathname.
 * @author Fischl
 */
public final class ClusteringShardWriter{

	private ClusteringShardWriter(){}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	/** Preserves original headers/residue bytes, with appendFasta line wrapping.
	 * All outputs must have existing writable parent directories. No short-record
	 * filtering occurs. Empty input produces empty shards and a mapping header.
	 */
	public static void write(ClusteringInput input, String fastaPattern, String mappingPath,
			int ways, ClusteringOrder.Mode orderMode, long seed, boolean overwrite){
		write(input, fastaPattern, mappingPath, ways, orderMode, seed, overwrite, null);
	}

	/** Optionally publishes one immutable receipt after each complete shard closes.
	 * A receipt permits reading that FASTA while later shards are written; it does
	 * not certify the still-open mapping or the whole partition. Receipt paths
	 * must be new even with overwrite enabled, preventing stale readiness reuse.
	 */
	public static void write(ClusteringInput input, String fastaPattern, String mappingPath,
			int ways, ClusteringOrder.Mode orderMode, long seed, boolean overwrite, String readyPattern){
		if(input==null || ways<1 ||
			(orderMode!=ClusteringOrder.Mode.SHUFFLE && orderMode!=ClusteringOrder.Mode.HYDRO)){
			throw new IllegalArgumentException("Sharding needs input, positive ways, and SHUFFLE or HYDRO order");
		}
		if(fastaPattern==null || fastaPattern.indexOf('#')<0 ||
			fastaPattern.indexOf('#')!=fastaPattern.lastIndexOf('#')){
			throw new IllegalArgumentException("Shard FASTA pattern must contain exactly one #");
		}
		if(readyPattern!=null && (readyPattern.indexOf('#')<0 ||
			readyPattern.indexOf('#')!=readyPattern.lastIndexOf('#'))){
			throw new IllegalArgumentException("Shard ready pattern must contain exactly one #");
		}
		final Path[] outputs=outputs(fastaPattern, mappingPath, ways, overwrite, readyPattern);
		final int[] order=permutation(input.sequences, orderMode, seed);
		final ByteBuilder fastaRow=new ByteBuilder(), mappingRow=new ByteBuilder();
		// The loader keeps the preservation flag private. Probe before opening any
		// output; an empty input has no original records whose bytes could be lost.
		if(order.length>0){input.appendFasta(0, fastaRow); fastaRow.clear();}
		try(Output mapping=new Output(outputs[ways], overwrite)){
			mapping.writer.println(MAPPING_HEADER);
			int emitted=0;
			for(int shard=0; shard<ways; shard++){
				final int perShard=shardSize(order.length, ways, shard);
				try(Output fasta=new Output(outputs[shard], overwrite)){
					for(int local=0; local<perShard; local++){
						final int parent=order[emitted++];
						fastaRow.clear(); input.appendFasta(parent, fastaRow); fasta.writer.print(fastaRow);
						mappingRow.clear().append(parent).tab().append(shard).tab().append(local).tab()
							.append(input.sequences[parent].id).nl();
						mapping.writer.print(mappingRow);
					}
				}
				// Deliberately outside close/finally: those also run on a partial
				// producer failure, which must never publish a readiness receipt.
				if(readyPattern!=null){
					publishReady(outputs[ways+1+shard], outputs[shard], shard, ways, perShard, order.length);
				}
			}
			if(emitted!=order.length){throw new IllegalStateException("Shard blocks must emit every validated parent ordinal exactly once");}
		}
	}

	/** Publishes complete metadata atomically without replacing an existing receipt.
	 * A same-directory hard link makes the closed temporary file visible in one
	 * operation and fails if the target exists; atomic move can replace a target.
	 */
	private static void publishReady(Path receipt, Path fasta, int shard, int ways, int records, int population){
		assert(records==shardSize(population, ways, shard)) :
			"Ready receipts must bind the complete expected shard, not a partially emitted prefix";
		Path temporary=null;
		try{
			temporary=Files.createTempFile(receipt.getParent(), ".shard-ready-", ".tmp");
			final String row="SHARD_READY_V1\t"+shard+"\t"+ways+"\t"+records+"\t"+population+"\t"+fasta+"\n";
			Files.write(temporary, row.getBytes(StandardCharsets.UTF_8));
			Files.createLink(receipt, temporary);
		}catch(IOException problem){throw new IllegalStateException("Could not publish completed shard receipt: "+receipt, problem);}
		finally{
			if(temporary!=null){
				try{Files.delete(temporary);}
				catch(IOException problem){throw new IllegalStateException("Could not remove shard receipt temporary: "+temporary, problem);}
			}
		}
	}

	/** Assigns the remainder to the first shards; divisible inputs retain the old boundaries. */
	public static int shardSize(int records, int ways, int shard){
		if(records<0 || ways<1 || shard<0 || shard>=ways){
			throw new IllegalArgumentException("Shard size needs nonnegative records, positive ways and an in-range shard");
		}
		final int size=records/ways+(shard<records%ways ? 1 : 0);
		assert(size>=0 && size<=records) : "Quotient/remainder blocks cannot contain more than the complete input";
		return size;
	}

	/** Matches the clusterer shuffle, then uses the existing stable hydrophobicity sort. */
	private static int[] permutation(ProteinSequence[] sequences, ClusteringOrder.Mode mode, long seed){
		final int[] order=new int[sequences.length];
		for(int i=0; i<order.length; i++){
			final ProteinSequence sequence=sequences[i];
			if(sequence==null || sequence.id==null || sequence.id.isEmpty() || sequence.id.indexOf('\t')>=0 ||
				sequence.id.indexOf('\n')>=0 || sequence.id.indexOf('\r')>=0){
				throw new IllegalArgumentException("Missing input or unsafe mapping identifier at parent ordinal "+i);
			}
			order[i]=i;
		}
		final Random random=new Random(seed);
		for(int i=order.length-1; i>0; i--){
			final int j=random.nextInt(i+1), value=order[i]; order[i]=order[j]; order[j]=value;
		}
		ClusteringOrder.apply(order, sequences, mode);
		final BitSet seen=new BitSet(order.length);
		for(int parent:order){
			if(parent<0 || parent>=order.length || seen.get(parent)){
				throw new IllegalStateException("Shard ordering must preserve disjoint parent-ordinal coverage");
			}
			seen.set(parent);
		}
		assert(seen.cardinality()==order.length) : "A distinct in-range permutation of all input ordinals must cover the complete input";
		return order;
	}

	/** Resolves every output before any writer opens, including symlink and hard-link aliases. */
	private static Path[] outputs(String pattern, String mapping, int ways, boolean overwrite, String readyPattern){
		final Path[] paths=new Path[Math.addExact(ways, readyPattern==null ? 1 : Math.addExact(ways, 1))];
		final HashSet<Path> canonical=new HashSet<>();
		final ArrayList<Path> existing=new ArrayList<>();
		try{
			for(int i=0; i<paths.length; i++){
				final boolean receipt=i>ways;
				final String name=i==ways ? mapping : receipt ?
					readyPattern.replace("#", Integer.toString(i-ways-1)) : pattern.replace("#", Integer.toString(i));
				if(name==null || name.isEmpty()){throw new IllegalArgumentException("Every shard and mapping needs a disk filename");}
				if(readyPattern!=null && (name.indexOf('\t')>=0 || name.indexOf('\n')>=0 || name.indexOf('\r')>=0)){
					throw new IllegalArgumentException("Ready receipt paths must fit one tab-delimited row: "+name);
				}
				final String lower=name.toLowerCase(Locale.ROOT);
				for(String stream:new String[]{"stdin", "stdout", "stderr"}){
					if(lower.equals(stream) || lower.startsWith(stream+".")){
						throw new IllegalArgumentException("Shards and mapping require disk files, not standard streams: "+name);
					}
				}
				final Path requested=Paths.get(name).toAbsolutePath().normalize();
				final Path parent=requested.getParent();
				if(parent==null || !Files.isDirectory(parent) || !Files.isWritable(parent)){
					throw new IllegalArgumentException("Output parent must already be a writable directory: "+name);
				}
				final boolean present=Files.exists(requested, LinkOption.NOFOLLOW_LINKS);
				if(present && (receipt || !overwrite || !Files.isRegularFile(requested) || !Files.isWritable(requested))){
					throw new IllegalArgumentException("Output exists without overwrite or is not a writable regular file: "+name);
				}
				final Path path=present ? requested.toRealPath() : parent.toRealPath().resolve(requested.getFileName());
				if(!canonical.add(path)){throw new IllegalArgumentException("Shard outputs alias the same path: "+name);}
				if(present){
					// Only pre-existing overwrite targets need inode comparisons; fresh
					// output sets use the linear canonical-path set above.
					for(Path previous:existing){
						if(Files.isSameFile(previous, path)){throw new IllegalArgumentException("Shard outputs alias the same file: "+name);}
					}
					existing.add(path);
				}
				// Resolve identity for collision checks, but preserve the requested
				// filename suffix for native format/compression detection through links.
				paths[i]=requested;
			}
		}catch(IOException problem){throw new IllegalArgumentException("Could not resolve shard output paths", problem);}
		return paths;
	}

	/** Native writer lifetime, joined even when a sibling writer or producer fails. */
	private static final class Output implements AutoCloseable{
		Output(Path path, boolean overwrite){
			this.path=path;
			writer=new ByteStreamWriter(path.toString(), overwrite, false, true);
			writer.start();
		}
		@Override
		public void close(){
			if(writer.poisonAndWait()){throw new IllegalStateException("Failed shard output: "+path);}
		}
		final Path path;
		final ByteStreamWriter writer;
	}

	public static final String MAPPING_HEADER="parent_input_index\tshard\tshard_input_index\tinput_id";
}
