package assemble;

import java.io.File;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectIntMap;
import parse.LineParser1;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Streams same-K features and their separate candidate/control truth namespaces
 * into native trainer tables. Unknown truth is counted, never made negative.
 * Genome folds are supplied by a frozen manifest, never selected from labels.
 * @author Fischl
 */
public final class SameKTrainingCorpus {

	/** Runs the converter; success requires a complete marker after every writer closes. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, SameKTrainingCorpus.class, false);
		try{new SameKTrainingCorpus(pp.args, pp.outstream).run();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Requires frozen cohort/inventory inputs, a collection root and a new output directory. */
	private SameKTrainingCorpus(final String[] args, final PrintStream stream){
		outstream=stream;
		final Parser parser=new Parser();
		for(String arg : args){
			final String[] split=arg.split("=", 2);
			String name=split[0].toLowerCase(Locale.ROOT);
			while(name.startsWith("-")){name=name.substring(1);}
			final String value=split.length>1 ? split[1] : null;
			if(name.equals("root")){root=value;}
			else if(name.equals("inventory")){inventory=value;}
			else if(name.equals("expected")){expected=Integer.parseInt(value);}
			else if(!parser.parse(arg, name, value)){throw new IllegalArgumentException("Unknown corpus option: "+arg);}
		}
		input=parser.in1; output=parser.out1;
		if(input==null || output==null || root==null || inventory==null || expected<1 ||
				parser.in2!=null || parser.out2!=null){
			throw new IllegalArgumentException("Required in=cohort.tsv inventory=inventory.tsv root=collection_root out=new_directory.");
		}
		if(new File(output).exists() || !Tools.testInputFiles(false, true, input, inventory)){
			throw new IllegalArgumentException("Inputs must exist and output directory must be new.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Preflights identities/files, streams each genome once, then closes all outputs. */
	private void run(){
		final ArrayList<Sample> samples=readSamples();
		checkInventory(samples);
		for(Sample sample : samples){
			final String dir=root+"/runs/"+sample.id;
			for(String file : FILES){
				if(!Tools.testInputFiles(false, true, dir+"/"+file)){
					throw new IllegalArgumentException("Missing source file: "+dir+"/"+file);
				}
			}
			final ByteFile exit=ByteFile.makeByteFile(dir+"/exit_code.txt", true);
			try{require(exit.nextLine(), "0"); if(exit.nextLine()!=null){throw new IllegalArgumentException("Extra exit-code data");}}
			finally{if(exit.close()){throw new IllegalStateException("Could not read exit code: "+dir);}}
		}
		if(!new File(output).mkdir()){throw new IllegalArgumentException("Cannot create new output directory: "+output);}
		final ByteStreamWriter[] writers=new ByteStreamWriter[5];
		final long[] emitted=new long[3], totals=new long[6];
		boolean error=false;
		try{
			for(int i=0; i<3; i++){
				writers[i]=open(output+"/"+FOLDS[i]+".tsv");
				writers[i].print("#dims\t27\t1\n#schema=samek27-v1; unknown labels excluded; no class balancing\n");
			}
			writers[3]=open(output+"/rows.tsv");
			writers[3].print("fold\trow\tgenome\tkind\tevent\ttarget\n");
			writers[4]=open(output+"/census.tsv");
			writers[4].print("genome\tgroup\tfold\tdepth\tk\tcandidate_supported\tcandidate_contradicted\tcandidate_na\tcontrol_supported\tcontrol_contradicted\tcontrol_na\tcontrol_population\tcontrol_sampled\tcontrol_requested\tcontrol_seed\treference_sha80\n");
			for(Sample sample : samples){
				final long[] counts=convert(sample, writers, emitted);
				for(int i=0; i<6; i++){totals[i]+=counts[i];}
			}
		}finally{
			for(ByteStreamWriter writer : writers){if(writer!=null){error=writer.poisonAndWait() | error;}}
		}
		if(error){throw new IllegalStateException("Incomplete training-table output.");}
		final long resolved=totals[0]+totals[1]+totals[3]+totals[4];
		assert(resolved==emitted[0]+emitted[1]+emitted[2]) : "Each known target emits exactly one trainer row and provenance row.";
		final String summary="SAME_K_TRAINING_CORPUS_PASS genomes="+samples.size()+" emitted="+resolved+
				" train="+emitted[0]+" validation="+emitted[1]+" test="+emitted[2]+
				" supported="+(totals[0]+totals[3])+" contradicted="+(totals[1]+totals[4])+" unresolved="+(totals[2]+totals[5]);
		final ByteStreamWriter complete=open(output+"/complete.txt");
		complete.print(summary+"\n");
		if(complete.poisonAndWait()){throw new IllegalStateException("Could not save completion marker.");}
		outstream.println(summary);
	}

	/** Validates unique genomes/taxa/paths, selected K membership and group-disjoint folds. */
	private ArrayList<Sample> readSamples(){
		assert(expected>0) : "A caller-specified expected genome count detects truncated manifests.";
		final ArrayList<Sample> samples=new ArrayList<Sample>();
		final ObjectIntMap<String> ids=new ObjectIntMap<String>(String.class), groups=new ObjectIntMap<String>(String.class);
		final ObjectIntMap<String> taxa=new ObjectIntMap<String>(String.class), paths=new ObjectIntMap<String>(String.class);
		final ByteFile file=ByteFile.makeByteFile(input, true);
		final LineParser1 fields=new LineParser1('\t');
		try{
			final byte[] metadata=file.nextLine();
			if(metadata==null || !new String(metadata, StandardCharsets.UTF_8).startsWith("#samek_corpus_v1 ")){
				throw new IllegalArgumentException("Missing frozen cohort version.");
			}
			require(file.nextLine(), COHORT_HEADER);
			for(byte[] line=file.nextLine(); line!=null; line=file.nextLine()){
				fields.set(line);
				if(fields.terms()!=16){throw new IllegalArgumentException("Expected sixteen cohort fields.");}
				final String id=fields.parseString(0), group=fields.parseString(3), taxid=fields.parseString(2);
				final String path=fields.parseString(1), fold=fields.parseString(12);
				int split=-1;
				for(int i=0; i<3; i++){if(FOLDS[i].equals(fold)){split=i;}}
				final int k=fields.parseInt(13), depth=fields.parseInt(7);
				final long seed=fields.parseLong(10);
				if(!id.matches("g[0-9]{4}") || group.isEmpty() || split<0 || k<24 || k>127 ||
						depth<15 || depth>150 || seed<0 || ids.get(id)>=0 || taxa.get(taxid)>=0 || paths.get(path)>=0 ||
						(groups.get(group)>=0 && groups.get(group)!=split)){
					throw new IllegalArgumentException("Invalid/duplicate identity or leaking group at "+id);
				}
				boolean found=false;
				for(String token : fields.parseString(9).split(",")){found|=Integer.parseInt(token)==k;}
				if(!found){throw new IllegalArgumentException("K is not in original schedule: "+id);}
				ids.put(id, split); taxa.put(taxid, split); paths.put(path, split); groups.put(group, split);
				samples.add(new Sample(id, group, split, k, depth, seed));
			}
		}finally{if(file.close()){throw new IllegalStateException("Cohort read failure.");}}
		if(samples.size()!=expected){throw new IllegalArgumentException("Genome count "+samples.size()+" != "+expected);}
		return samples;
	}

	/** Binds each cohort member to the prior reference check and rejects byte-identical references. */
	private void checkInventory(final ArrayList<Sample> samples){
		assert(samples.size()==expected) : "Inventory is compared to the complete frozen cohort in original order.";
		final ByteFile file=ByteFile.makeByteFile(inventory, true);
		final LineParser1 fields=new LineParser1('\t');
		final ObjectIntMap<String> hashes=new ObjectIntMap<String>(String.class);
		try{
			require(file.nextLine(), "genome\tread_bytes\treference_bytes\tverified_reference_sha80");
			for(Sample sample : samples){
				final byte[] line=file.nextLine();
				if(line==null){throw new IllegalArgumentException("Truncated inventory at "+sample.id);}
				fields.set(line);
				if(fields.terms()!=4 || !fields.termEquals(sample.id, 0) || fields.parseLong(1)<1 || fields.parseLong(2)<1){
					throw new IllegalArgumentException("Inventory identity/size mismatch at "+sample.id);
				}
				sample.hash=fields.parseString(3);
				if(!sample.hash.matches("[0-9a-f]{20}") || hashes.get(sample.hash)>=0){
					throw new IllegalArgumentException("Invalid or duplicate reference hash at "+sample.id);
				}
				hashes.put(sample.hash, sample.fold);
			}
			if(file.nextLine()!=null){throw new IllegalArgumentException("Extra inventory rows.");}
		}finally{if(file.close()){throw new IllegalStateException("Inventory read failure.");}}
	}

	/** Joins all feature/context/truth rows in lockstep with O(1) per-row scratch. */
	private long[] convert(final Sample sample, final ByteStreamWriter[] out, final long[] emitted){
		assert(sample.fold>=0 && sample.fold<3) : "Frozen genome-level folds route rows; labels never select a fold.";
		final String dir=root+"/runs/"+sample.id;
		final ByteFile[] files=new ByteFile[5];
		final long[] seen=new long[2], counts=new long[6];
		final LineParser1 fp=new LineParser1('\t'), cp=new LineParser1('\t'), lp=new LineParser1('\t');
		final ByteBuilder row=new ByteBuilder();
		boolean trailer=false, inControls=false, error=false;
		long population=-1, sampled=-1, requested=-1;
		try{
			for(int i=0; i<5; i++){files[i]=ByteFile.makeByteFile(dir+"/"+FILES[i], true);}
			final String[] meta=text(files[0].nextLine()).split(" ", -1);
			if(meta.length!=4 || !meta[0].equals("#schema="+SameKJoinFeatures.VERSION)){throw new IllegalArgumentException("Wrong feature schema.");}
			final long unique=number(meta[1], "unique=");
			final double observations=decimal(meta[2], "observations="), squares=decimal(meta[3], "squares=");
			if(unique<0 || observations<unique || squares<observations || (unique==0 && (observations!=0 || squares!=0))){
				throw new IllegalArgumentException("Invalid depth moment metadata.");
			}
			final float mean=SameKJoinFeatures.depth(unique==0 ? 0 : observations/unique);
			final float weighted=SameKJoinFeatures.depth(observations==0 ? 0 : squares/observations);
			final float encodedK=(float)(0.125*Math.log(sample.k)*(1/Math.log(2)));
			final ByteBuilder header=new ByteBuilder("kind\tevent");
			for(int i=0; i<27; i++){header.tab().append('f').append(i);}
			require(files[0].nextLine(), header.toString());
			require(files[1].nextLine(), SameKJoinDiagnostic.HEADER);
			require(files[2].nextLine(), SameKJoinExport.CONTROL_HEADER);
			for(int i=0; i<2; i++){
				require(files[3+i].nextLine(), LABEL_COMMENT);
				require(files[3+i].nextLine(), (i==0 ? "event" : "control_event")+LABEL_SUFFIX);
			}
			for(byte[] line=files[0].nextLine(); line!=null; line=files[0].nextLine()){
				if(line.length>0 && line[0]=='#'){
					if(trailer){throw new IllegalArgumentException("Duplicate sampling trailer.");}
					final String[] terms=text(line).split(" ", -1);
					if(terms.length!=5 || !terms[0].equals("#controls")){throw new IllegalArgumentException("Unexpected feature comment.");}
					population=number(terms[1], "population="); sampled=number(terms[2], "sampled=");
					if(number(terms[3], "seed=")!=sample.seed){throw new IllegalArgumentException("Control seed mismatch at "+sample.id);}
					requested=number(terms[4], "requested="); trailer=true;
					continue;
				}
				if(trailer){throw new IllegalArgumentException("Feature row after sampling trailer.");}
				fp.set(line);
				if(fp.terms()!=29){throw new IllegalArgumentException("Expected kind/event plus27 features at "+sample.id);}
				final int kind=fp.termEquals("candidate", 0) ? 0 : fp.termEquals("control", 0) ? 1 : -1;
				if(kind<0 || (kind==0 && inControls)){throw new IllegalArgumentException("Invalid population order.");}
				inControls|=kind==1;
				final long event=++seen[kind];
				if(fp.parseLong(1)!=event){throw new IllegalArgumentException("Nonconsecutive feature ID at "+sample.id);}
				final byte[] context=files[1+kind].nextLine(), label=files[3+kind].nextLine();
				if(context==null || label==null){throw new IllegalArgumentException("Missing paired context/label at "+sample.id);}
				cp.set(context); lp.set(label);
				if(cp.terms()!=(kind==0 ? 36 : 10) || lp.terms()!=17 || cp.parseLong(0)!=event || lp.parseLong(0)!=event ||
						cp.parseInt(1)!=sample.k || lp.parseInt(1)!=sample.k){
					throw new IllegalArgumentException("Context/label identity or K mismatch at "+sample.id+" event="+event);
				}
				for(int i=2; i<29; i++){
					final float value=fp.parseFloat(i);
					if(!Float.isFinite(value) || value<0){throw new IllegalArgumentException("Invalid feature at "+sample.id);}
				}
				if(fp.parseFloat(26)!=encodedK || fp.parseFloat(27)!=mean || fp.parseFloat(28)!=weighted){
					throw new IllegalArgumentException("Feature K/global moments disagree with provenance at "+sample.id);
				}
				final int labelClass=lp.termEquals("supported", 2) && lp.termEquals("1", 3) ? 0 :
						lp.termEquals("contradicted", 2) && lp.termEquals("0", 3) ? 1 :
						lp.termEquals("unresolved", 2) && lp.termEquals("NA", 3) ? 2 : -1;
				if(labelClass<0){throw new IllegalArgumentException("Invalid truth/target pair at "+sample.id);}
				counts[kind*3+labelClass]++;
				if(labelClass==2){continue;}
				int start=0;
				for(int tabs=0; tabs<2; start++){if(line[start]=='\t'){tabs++;}}
				final int target=1-labelClass;
				out[sample.fold].print(row.clear().append(line, start, line.length-start).tab().append(target).nl());
				out[3].print(row.clear().append(FOLDS[sample.fold]).tab().append(++emitted[sample.fold]).tab().append(sample.id)
						.tab().append(kind==0 ? "candidate" : "control").tab().append(event).tab().append(target).nl());
			}
			for(int i=1; i<5; i++){if(files[i].nextLine()!=null){throw new IllegalArgumentException("Unpaired extra source rows at "+sample.id);}}
			if(!trailer || population<0 || requested<0 || sampled!=seen[1] || sampled!=Math.min(population, requested)){
				throw new IllegalArgumentException("Control sampling/count mismatch at "+sample.id);
			}
			row.clear().append(sample.id).tab().append(sample.group).tab().append(FOLDS[sample.fold]).tab().append(sample.depth).tab().append(sample.k);
			for(long count : counts){row.tab().append(count);}
			out[4].print(row.tab().append(population).tab().append(sampled).tab().append(requested).tab().append(sample.seed).tab().append(sample.hash).nl());
		}finally{
			for(ByteFile file : files){if(file!=null){error=file.close() | error;}}
		}
		if(error){throw new IllegalStateException("Incomplete corpus input at "+sample.id);}
		return counts;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Opens a fresh output writer; callers always close it even after validation failure. */
	private static ByteStreamWriter open(final String path){
		if(!Tools.testOutputFiles(false, false, false, path)){throw new IllegalArgumentException("Output exists: "+path);}
		final ByteStreamWriter writer=new ByteStreamWriter(path, false, false, false);
		writer.start();
		return writer;
	}

	/** Exact headers keep similarly shaped future schemas from silently entering training. */
	private static void require(final byte[] line, final String expected){
		if(!Arrays.equals(line, expected.getBytes(StandardCharsets.UTF_8))){throw new IllegalArgumentException("Unexpected corpus header; expected "+expected);}
	}

	/** Only cold metadata becomes a String; all feature rows stay byte-oriented. */
	private static String text(final byte[] line){
		if(line==null){throw new IllegalArgumentException("Missing metadata line.");}
		return new String(line, StandardCharsets.UTF_8);
	}

	/** Checks a named integral metadata field, including a full signed-long parse. */
	private static long number(final String field, final String prefix){
		if(!field.startsWith(prefix)){throw new IllegalArgumentException("Expected metadata field "+prefix);}
		return Long.parseLong(field.substring(prefix.length()));
	}

	/** Metadata doubles may use exponents; native feature decimals never do. */
	private static double decimal(final String field, final String prefix){
		if(!field.startsWith(prefix)){throw new IllegalArgumentException("Expected metadata field "+prefix);}
		final double value=Double.parseDouble(field.substring(prefix.length()));
		if(!Double.isFinite(value) || value<0){throw new IllegalArgumentException("Nonfinite/negative moment.");}
		return value;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Retains one genome's small routing record, not its feature vectors. */
	private static final class Sample {
		/** Stores the label-independent route and the expected K/control seed. */
		Sample(final String id_, final String group_, final int fold_, final int k_, final int depth_, final long seed_){
			id=id_; group=group_; fold=fold_; k=k_; depth=depth_; seed=seed_;
		}
		final String id, group;
		final int fold, k, depth;
		final long seed;
		String hash;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final PrintStream outstream;
	private String input, inventory, root, output;
	private int expected=1000;
	private static final String[] FOLDS={"train", "validation", "test"};
	private static final String[] FILES={"features.tsv", "candidates.tsv", "controls.tsv", "labels.tsv", "control_labels.tsv", "complete.txt", "exit_code.txt"};
	private static final String COHORT_HEADER="#genome\tsource\ttaxid\tgroup\tcontig_length50\tbases\tscaffolds\tdepth\tquality\tk_schedule\tread_seed\tplan_seed\tfold\tk\treads\treference";
	private static final String LABEL_COMMENT="#pair labels use bounded exact context, not whole-assembly truth; ambiguous placements remain unresolved";
	private static final String LABEL_SUFFIX="\tk\tlabel\ttarget\twhole_hits\tleft_hits\tright_hits\tsource_hits\tsource_unknown\tsource_left_mask\tsource_right_mask\tsource_pair_mask\tdest_hits\tdest_unknown\tdest_left_mask\tdest_right_mask\tdest_pair_mask";
}
