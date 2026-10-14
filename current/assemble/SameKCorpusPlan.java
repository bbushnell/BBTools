package assemble;

import java.io.File;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Locale;
import java.util.Random;

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
 * Freezes single-K collection recipes before labels are examined. All versions
 * of a recorded source group share one fold; groups are filename-derived, not
 * verified genera. This planner never reads reference sequence or training labels.
 * @author Fischl
 */
public final class SameKCorpusPlan {

	/** Parses a native BBTools command and writes a new collection manifest. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, SameKCorpusPlan.class, false);
		try{new SameKCorpusPlan(pp.args, pp.outstream).run();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Accepts a frozen Plan35 sample manifest and the existing library root. */
	private SameKCorpusPlan(final String[] args, final PrintStream stream){
		outstream=stream;
		final Parser parser=new Parser();
		for(String arg : args){
			final String[] split=arg.split("=", 2);
			String name=split[0].toLowerCase(Locale.ROOT);
			while(name.startsWith("-")){name=name.substring(1);}
			final String value=split.length>1 ? split[1] : null;
			if(name.equals("root")){root=value;}
			else if(name.equals("splitseed")){splitSeed=Long.parseLong(value);}
			else if(name.equals("kseed")){kSeed=Long.parseLong(value);}
			else if(name.equals("expected")){expected=Integer.parseInt(value);}
			else if(!parser.parse(arg, name, value)){
				throw new IllegalArgumentException("Unknown corpus-plan argument: "+arg);
			}
		}
		input=parser.in1;
		output=parser.out1;
		if(input==null || output==null || root==null || expected<3 ||
				parser.in2!=null || parser.out2!=null || !new File(root).isAbsolute()){
			throw new IllegalArgumentException("Required in=samples.tsv out=new.tsv root=/absolute/library/root expected=N>=3.");
		}
		TadpoleGraph.checkPaths(input, new ArrayList<String>(), output);
		if(!Tools.testInputFiles(false, true, input) || !Tools.testOutputFiles(false, false, false, output)){
			throw new IllegalArgumentException("Input must exist; manifest output must be new.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Validates everything before opening output, then freezes folds and K draws. */
	private void run(){
		assert(expected>=3) : "Three nonempty folds require at least three source groups.";
		final ArrayList<Sample> samples=readSamples();
		final ObjectIntMap<String> folds=new ObjectIntMap<String>(String.class);
		final ArrayList<String> groups=new ArrayList<String>();
		for(Sample sample : samples){
			if(folds.get(sample.group)<0){
				folds.put(sample.group, 0);
				groups.add(sample.group);
			}
		}
		if(groups.size()<3){throw new IllegalArgumentException("Need at least three distinct groups for disjoint folds.");}
		// Sort before seeded shuffling so the split does not depend on input order.
		Collections.sort(groups);
		Collections.shuffle(groups, new Random(splitSeed));
		final int heldout=Math.max(1, (int)Math.round(groups.size()*0.15));
		final int trainEnd=groups.size()-2*heldout, validEnd=groups.size()-heldout;
		for(int i=0; i<groups.size(); i++){
			folds.put(groups.get(i), i<trainEnd ? 0 : i<validEnd ? 1 : 2);
		}
		// Stable genome ordering also makes each random K draw independent of row order.
		Collections.sort(samples, (a, b) -> a.id.compareTo(b.id));
		final Random random=new Random(kSeed);
		final int[] rows=new int[3];
		for(Sample sample : samples){
			sample.fold=folds.get(sample.group);
			sample.k=sample.ks[random.nextInt(sample.ks.length)];
			rows[sample.fold]++;
		}
		write(samples, groups.size(), trainEnd, heldout);
		outstream.println("SAME_K_CORPUS_PLAN_PASS genomes="+samples.size()+" groups="+groups.size()+
				" train="+rows[0]+" validation="+rows[1]+" test="+rows[2]+
				" splitseed="+splitSeed+" kseed="+kSeed);
	}

	/** Reads only the small manifest; validates source identity and simulation ranges. */
	private ArrayList<Sample> readSamples(){
		assert(input!=null) : "A frozen source manifest is required for reproducible library reuse.";
		final ArrayList<Sample> samples=new ArrayList<Sample>();
		final ObjectIntMap<String> ids=new ObjectIntMap<String>(String.class);
		final ObjectIntMap<String> paths=new ObjectIntMap<String>(String.class);
		final ObjectIntMap<String> taxa=new ObjectIntMap<String>(String.class);
		final ByteFile file=ByteFile.makeByteFile(input, true);
		final LineParser1 fields=new LineParser1('\t');
		boolean error;
		try{
			if(!Arrays.equals(file.nextLine(), HEADER.getBytes(StandardCharsets.UTF_8))){
				throw new IllegalArgumentException("Expected the exact twelve-column Plan35 samples header.");
			}
			for(byte[] line=file.nextLine(); line!=null; line=file.nextLine()){
				fields.set(line);
				if(fields.terms()!=12){throw new IllegalArgumentException("Expected twelve fields at row "+(samples.size()+1));}
				final String id=fields.parseString(0), path=fields.parseString(1);
				final String taxid=fields.parseString(2), group=fields.parseString(3);
				if(!id.matches("g[0-9]{4}") || !new File(path).isAbsolute() || group.length()==0 ||
						ids.get(id)>=0 || paths.get(path)>=0 || taxa.get(taxid)>=0){
					throw new IllegalArgumentException("Duplicate or invalid genome/path/taxid/group at "+id);
				}
				if(fields.parseLong(2)<1 || fields.parseLong(4)<1 || fields.parseLong(4)>fields.parseLong(5) ||
						fields.parseLong(5)>15000000 || fields.parseInt(6)<1 ||
						fields.parseInt(7)<15 || fields.parseInt(7)>150 ||
						fields.parseInt(8)<25 || fields.parseInt(8)>35 ||
						fields.parseLong(10)<0 || fields.parseLong(11)!=350194){
					throw new IllegalArgumentException("Unexpected original simulation recipe at "+id);
				}
				final String[] tokens=fields.parseString(9).split(",", -1);
				if(tokens.length<3 || tokens.length>4){throw new IllegalArgumentException("Expected three or four original Ks at "+id);}
				final int[] ks=new int[tokens.length];
				for(int i=0; i<ks.length; i++){
					ks[i]=Integer.parseInt(tokens[i]);
					if(ks[i]<24 || ks[i]>127 || (i>0 && ks[i]<=ks[i-1])){
						throw new IllegalArgumentException("Ks must be distinct, ascending and in [24,127] at "+id);
					}
				}
				ids.put(id, samples.size()); paths.put(path, samples.size()); taxa.put(taxid, samples.size());
				samples.add(new Sample(line, id, group, ks));
			}
		}finally{error=file.close();}
		if(error || samples.size()!=expected){
			throw new IllegalArgumentException("Incomplete source manifest: rows="+samples.size()+", expected="+expected);
		}
		return samples;
	}

	/** Writes original provenance unchanged, followed by fold, K and exact library paths. */
	private void write(final ArrayList<Sample> samples, final int groups, final int trainGroups, final int heldout){
		assert(samples.size()==expected) : "Only a fully validated manifest may be emitted.";
		final ByteStreamWriter writer=new ByteStreamWriter(output, false, false, false);
		writer.start();
		boolean error;
		try{
			writer.print("#samek_corpus_v1 splitseed="+splitSeed+" kseed="+kSeed+
					" groups="+groups+" train_groups="+trainGroups+" validation_groups="+heldout+
					" test_groups="+heldout+" group_source=filename input_files_verified=false\n");
			writer.print(HEADER+"\tfold\tk\treads\treference\n");
			final ByteBuilder row=new ByteBuilder();
			for(Sample sample : samples){
				final String directory=new File(root, sample.id).getPath();
				writer.print(row.clear().append(sample.original).tab().append(FOLDS[sample.fold]).tab().append(sample.k)
						.tab().append(directory).append("/reads.fq.gz").tab().append(directory).append("/reference.fa").nl());
			}
		}finally{error=writer.poisonAndWait();}
		if(error){throw new IllegalStateException("Failed to write complete corpus plan: "+output);}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Retains a small source row and its assigned collection recipe, never sequence. */
	private static final class Sample {
		Sample(final byte[] original_, final String id_, final String group_, final int[] ks_){
			assert(original_!=null && ks_.length>=3) : "Samples retain validated provenance and original K choices.";
			original=original_; id=id_; group=group_; ks=ks_;
		}
		final byte[] original;
		final String id, group;
		final int[] ks;
		int fold, k;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final PrintStream outstream;
	private String input, output, root;
	private int expected=1000;
	private long splitSeed=370194, kSeed=370195;
	private static final String[] FOLDS={"train", "validation", "test"};
	private static final String HEADER="#genome\tsource\ttaxid\tgroup\tcontig_length50\tbases\tscaffolds\tdepth\tquality\tk_schedule\tread_seed\tplan_seed";
}
