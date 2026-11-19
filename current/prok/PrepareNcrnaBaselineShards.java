package prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Locale;

import map.IntObjectMap;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import structures.ByteBuilder;

/** Partitions an accepted complete source-frame ledger without reading DNA.
 * Adapted from template.A_SampleBasic2's parser/output lifecycle and
 * PrepareNcrnaRosterBaseline's native byte prepass. Accepted assembly groups
 * remain together; largest groups are placed on the least-loaded shard using
 * their saved available nucleotide counts. Original frame names, rows, paths
 * and per-frame gene counts are unchanged. The full original-call ledger stays
 * the canonical grading population, including unavailable and zero-gene frames.
 * @author Brian Bushnell, Raiden
 */
public final class PrepareNcrnaBaselineShards{

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final PreParser pp=new PreParser(args, PrepareNcrnaBaselineShards.class, false);
		try{
			final PrepareNcrnaBaselineShards tool=new PrepareNcrnaBaselineShards(pp.args);
			tool.process();
			timer.stop();
			pp.outstream.println("NCRNA_BASELINE_SHARDS_COMPLETE time="+timer);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private PrepareNcrnaBaselineShards(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("shards")){shardCount=Parse.parseIntKMG(b);}
			else if(a.equals("expectedframes")){expectedFrames=Parse.parseIntKMG(b);}
			else if(a.equals("expectedgenes")){expectedGenes=Parse.parseIntKMG(b);}
			else{require(parser.parse(arg, a, b), "Unknown baseline sharding argument: "+arg);}
		}
		input=parser.in1;
		require(input!=null && parser.out1!=null && shardCount>0 && expectedFrames>=shardCount
			&& expectedGenes>=0 && !parser.overwrite && !parser.append && parser.in2==null
			&& parser.out2==null && parser.maxReads<0 && parser.samplerate==1,
			"Supply in=accepted.frames.tsv out=fresh_directory shards=N expectedframes=N expectedgenes=N without sampling");
		root=new File(parser.out1).getAbsoluteFile();
		require(!root.exists(), "Shard output must be fresh: "+root);
		Shared.setThreads(1);
	}

	private void process(){
		loadFrames();
		assignGroups();
		require(root.mkdir(), "Could not create fresh shard directory: "+root);
		for(final Shard shard:shards){writeShard(shard);}
		writePlan();
	}

	/*--------------------------------------------------------------*/
	/*----------------       Whole Population       ----------------*/
	/*--------------------------------------------------------------*/

	private void loadFrames(){
		final ObjectMap<String, Frame> names=new ObjectMap<String, Frame>(String.class, Frame.class);
		final ObjectMap<String, Frame> paths=new ObjectMap<String, Frame>(String.class, Frame.class);
		try(Input in=new Input(input)){
			RrnaResourceIO.header(in.reader, FRAME_HEADER, input);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				require(p.terms()==9, "Use the complete prepared nine-column frame ledger");
				final String name=p.parseString(0), path=p.parseString(2);
				final int groupId=p.parseInt(1), genes=p.parseInt(7);
				final long bases=p.parseLong(5);
				final boolean available=p.termEquals("AVAILABLE_RECORDED_FRAME", 8);
				require(available || p.termEquals("UNRESOLVED_SOURCE_UNAVAILABLE", 8), "Unknown source state cannot become a zero-call input");
				require(!name.isEmpty() && groupId>=0 && genes>=0 && bases>0 && p.parseLong(6)>0
					&& p.parseLong(4)>=0 && new File(path).isAbsolute() && path.indexOf(',')<0
					&& !names.contains(name) && !paths.contains(path), "Frame names and exact source paths must remain unique and valid");
				Group group=groupById.get(groupId);
				if(group==null){
					group=new Group(groupId);
					groupById.put(groupId, group);
					groups.add(group);
				}
				final Frame frame=new Frame(line, name, path, group, genes, available);
				names.put(name, frame);
				paths.put(path, frame);
				frames.add(frame);
				group.bases+=available ? bases : 0;
				totalGenes+=genes;
			}
		}
		require(frames.size()==expectedFrames && totalGenes==expectedGenes,
			"Shard planning must retain all recorded frames and their complete original-gene denominator");
	}

	/** Metadata-only largest-processing-time assignment; stable group/shard tie
	 * breaks make the emitted lists reproducible without rearranging frame rows. */
	private void assignGroups(){
		Collections.sort(groups, GROUP_ORDER);
		shards=new Shard[shardCount];
		for(int i=0; i<shards.length; i++){shards[i]=new Shard(i);}
		for(final Group group:groups){
			Shard best=shards[0];
			for(final Shard candidate:shards){
				if(candidate.bases<best.bases){best=candidate;}
			}
			group.shard=best;
			best.bases+=group.bases;
		}
		for(final Frame frame:frames){
			final Shard shard=frame.group.shard;
			shard.frames.add(frame);
			shard.genes+=frame.genes;
			if(frame.available){shard.inputs++;}
		}
		for(final Shard shard:shards){require(shard.inputs>0, "Each compute shard needs an available input; reduce shards= for this population");}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	private void writeShard(final Shard shard){
		assert(shard.inputs>0) : "CallGenes requires a nonempty explicit input list for each shard";
		final String prefix=prefix(shard);
		try(Output ledger=new Output(prefix+".frames.tsv"); Output reference=new Output(prefix+".reference.tsv");
			Output config=new Output(prefix+".inputs.config")){
			ledger.writer.println(FRAME_HEADER);
			reference.writer.println(GradeNcrnaReference.MANIFEST_HEADER);
			config.writer.print("in=");
			int written=0;
			for(final Frame frame:shard.frames){
				ledger.writer.print(row.clear().append(frame.raw).nl());
				if(!frame.available){continue;}
				reference.writer.print(row.clear().append(frame.name).tab().append(frame.path).tab()
					.append(prefix).append(".unused_").append(written).append(".gff\n"));
				row.clear();
				if(written>0){row.append(',');}
				config.writer.print(row.append(frame.path));
				written++;
			}
			config.writer.print("\n");
			require(written==shard.inputs, "Reference rows and caller input paths must agree including zero-gene inputs");
		}
	}

	/** Publish completion last, after every native output has closed. */
	private void writePlan(){
		long genes=0, count=0;
		try(Output ownership=new Output(new File(root, "ownership.tsv").getPath())){
			ownership.writer.println("shard_id\t"+FRAME_HEADER);
			for(final Frame frame:frames){
				ownership.writer.print(row.clear().append("shard_").append(frame.group.shard.id).tab().append(frame.raw).nl());
			}
		}
		try(Output plan=new Output(new File(root, "shards.tsv").getPath())){
			plan.writer.println("shard_id\trecorded_frames\tavailable_inputs\toriginal_genes\tavailable_bases\treference_manifest\tinput_config");
			for(final Shard shard:shards){
				count+=shard.frames.size();
				genes+=shard.genes;
				plan.writer.print(row.clear().append("shard_").append(shard.id).tab().append(shard.frames.size()).tab()
					.append(shard.inputs).tab().append(shard.genes).tab().append(shard.bases).tab()
					.append(prefix(shard)).append(".reference.tsv\t").append(prefix(shard)).append(".inputs.config\n"));
			}
			require(count==expectedFrames && genes==expectedGenes, "No frame or original gene may disappear during assignment");
		}
	}

	private String prefix(final Shard shard){return new File(root, "shard_"+shard.id).getPath();}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------             Types            ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Frame{
		Frame(byte[] raw_, String name_, String path_, Group group_, int genes_, boolean available_){
			raw=raw_;
			name=name_;
			path=path_;
			group=group_;
			genes=genes_;
			available=available_;
		}
		final byte[] raw;
		final String name, path;
		final Group group;
		final int genes;
		final boolean available;
	}
	private static final class Group{
		Group(int id_){id=id_;}
		final int id;
		long bases;
		Shard shard;
	}
	private static final class Shard{
		Shard(int id_){id=id_;}
		final int id;
		long bases, genes;
		int inputs;
		final ArrayList<Frame> frames=new ArrayList<Frame>();
	}
	private static final Comparator<Group> GROUP_ORDER=(a, b)->{
		final int order=Long.compare(b.bases, a.bases);
		return order==0 ? Integer.compare(a.id, b.id) : order;
	};
	private String input;
	private File root;
	private int shardCount, expectedFrames;
	private long expectedGenes=-1, totalGenes;
	private Shard[] shards;
	private final ArrayList<Frame> frames=new ArrayList<Frame>();
	private final ArrayList<Group> groups=new ArrayList<Group>();
	private final IntObjectMap<Group> groupById=new IntObjectMap<Group>(16, Group.class);
	private final LineParser1 p=new LineParser1('\t');
	private final ByteBuilder row=new ByteBuilder();
	private static final String FRAME_HEADER="input_id\tassembly_group\tcoordinate_frame\tsupplied_fasta\tfile_bytes\tnucleotide_bases\tfasta_records\toriginal_genes\tinput_state";
}
