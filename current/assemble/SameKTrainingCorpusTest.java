package assemble;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import structures.ByteBuilder;

/** Small hand-labeled schema fixtures; these test conversion, not biological truth. @author Fischl */
public final class SameKTrainingCorpusTest {

	/** Generates an isolated fixture or checks exact rows from a completed conversion. */
	public static void main(final String[] args){
		if(args.length<2){throw new IllegalArgumentException("generate root [variant], or check output");}
		if(args[0].equals("generate")){generate(args[1], args.length>2 ? args[2] : "valid");}
		else if(args[0].equals("check")){check(args[1]);}
		else{throw new IllegalArgumentException("Unknown fixture command.");}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates three genomes in three fixed folds, each with four candidates and two controls. */
	private static void generate(final String root, final String variant){
		if(!new File(root).mkdir() || !new File(root+"/runs").mkdir()){
			throw new IllegalArgumentException("Fixture directory must be new: "+root);
		}
		final ByteBuilder cohort=new ByteBuilder("#samek_corpus_v1 fixture=true\n");
		cohort.append("#genome\tsource\ttaxid\tgroup\tcontig_length50\tbases\tscaffolds\tdepth\tquality\tk_schedule\tread_seed\tplan_seed\tfold\tk\treads\treference\n");
		final ByteBuilder inventory=new ByteBuilder("genome\tread_bytes\treference_bytes\tverified_reference_sha80\n");
		for(int i=1; i<=3; i++){
			final String id="g000"+i, dir=root+"/runs/"+id;
			if(!new File(dir).mkdir()){throw new IllegalStateException("Cannot create fixture genome directory.");}
			final int k=32*i;
			final String mutation=i==1 ? variant : "valid";
			cohort.append(id).tab().append("/fixture/source").append(i).tab().append(i).tab().append("group")
					.append(variant.equals("group_leak") && i==2 ? 1 : i);
			cohort.append("\t1000\t1000\t1\t30\t30\t32,64,96\t").append(7000+i)
					.append("\t350194\t").append(FOLDS[i-1]).tab().append(k).append("\t/fixture/reads\t/fixture/reference\n");
			inventory.append(id).append("\t100\t100\t0000000000000000000")
					.append(variant.equals("duplicate_hash") && i==2 ? 1 : i).nl();
			final ByteBuilder features=new ByteBuilder("#schema=samek27-v1 unique=3 observations=6.0 squares=14.0\nkind\tevent");
			for(int j=0; j<27; j++){features.tab().append('f').append(j);}
			features.nl();
			final float[] vector=vector(k);
			for(int j=1; j<=4; j++){feature(features, "candidate", j, vector, j==1 ? mutation : "valid");}
			for(int j=1; j<=2; j++){feature(features, "control", j, vector, "valid");}
			if(!mutation.equals("missing_trailer")){
				features.append("#controls population=10 sampled=2 seed=")
						.append(mutation.equals("bad_seed") ? 0 : 7000+i).append(" requested=2\n");
			}
			write(dir+"/features.tsv", features);
			contexts(dir+"/candidates.tsv", SameKJoinDiagnostic.HEADER, 4, 36, k);
			contexts(dir+"/controls.tsv", SameKJoinExport.CONTROL_HEADER, 2, 10, k);
			labels(dir+"/labels.tsv", false, 4, k, mutation);
			labels(dir+"/control_labels.tsv", true, 2, k, mutation);
			write(dir+"/complete.txt", new ByteBuilder("fixture\n"));
			write(dir+"/exit_code.txt", new ByteBuilder(mutation.equals("failed_exit") ? "1\n" : "0\n"));
		}
		write(root+"/cohort.tsv", cohort);
		write(root+"/inventory.tsv", inventory);
	}

	/** Supplies finite normalized values through the already-tested shared feature producer. */
	private static float[] vector(final int k){
		final float[] values=new float[27];
		final int[] a={8, 7, 2, 1, 0, 6, 2, 0, 0}, b={9, 8, 3, 1, 0, 7, 3, 1, 0};
		SameKJoinFeatures.fill(values, k, a, b, 0, 0, 5, 6, 40, 41, 2, 14.0/6);
		return values;
	}

	/** Mutates exactly one boundary condition without corrupting unrelated fixture fields. */
	private static void feature(final ByteBuilder out, final String kind, final int event,
			final float[] values, final String mutation){
		out.append(kind).tab().append(mutation.equals("bad_id") ? event+1 : event);
		for(int i=0; i<values.length; i++){
			out.tab();
			if(i==0 && mutation.equals("nonfinite")){out.append("NaN");}
			else{out.append(i==25 && mutation.equals("bad_moment") ? values[i]+0.5f : values[i], 9);}
		}
		if(mutation.equals("bad_width")){out.append("\t1");}
		out.nl();
	}

	/** Context payload is deliberately dummy; this fixture checks schema and identity, not the oracle. */
	private static void contexts(final String path, final String header, final int count, final int width, final int k){
		final ByteBuilder out=new ByteBuilder(header).nl();
		for(int i=1; i<=count; i++){
			out.append(i).tab().append(k);
			for(int j=2; j<width; j++){out.append("\t0");}
			out.nl();
		}
		write(path, out);
	}

	/** Known labels cover both classes and unknowns; one mutation adds an unpaired label. */
	private static void labels(final String path, final boolean controls, final int count,
			final int k, final String mutation){
		final ByteBuilder out=new ByteBuilder("#pair labels use bounded exact context, not whole-assembly truth; ambiguous placements remain unresolved\n");
		out.append(controls && !mutation.equals("control_header") ? "control_event" : "event");
		out.append("\tk\tlabel\ttarget\twhole_hits\tleft_hits\tright_hits\tsource_hits\tsource_unknown\tsource_left_mask\tsource_right_mask\tsource_pair_mask\tdest_hits\tdest_unknown\tdest_left_mask\tdest_right_mask\tdest_pair_mask\n");
		final int n=count+(!controls && mutation.equals("extra_label") ? 1 : 0);
		for(int i=1; i<=n; i++){
			final boolean unknown=controls ? i==2 : i==3;
			final boolean negative=!controls && i==2;
			out.append(i).tab().append(i==1 && !controls && mutation.equals("bad_label_k") ? k+1 : k).tab();
			out.append(unknown ? "unresolved" : negative ? "contradicted" : "supported").tab();
			out.append(i==1 && !controls && mutation.equals("bad_target") ? "NA" : unknown ? "NA" : negative ? "0" : "1");
			for(int j=4; j<17; j++){out.append("\t0");}
			out.nl();
		}
		write(path, out);
	}

	/** Checks exact copied decimal rows, class counts and per-fold row provenance. */
	private static void check(final String root){
		final ByteBuilder expectedRows=new ByteBuilder("fold\trow\tgenome\tkind\tevent\ttarget\n");
		for(int i=0; i<3; i++){
			final ByteFile file=ByteFile.makeByteFile(root+"/"+FOLDS[i]+".tsv", true);
			try{
				equal(file.nextLine(), "#dims\t27\t1");
				if(file.nextLine()==null){throw new AssertionError("Missing schema comment");}
				for(int row=0; row<4; row++){
					final int target=row==1 ? 0 : 1;
					final ByteBuilder expected=new ByteBuilder();
					for(float value : vector(32*(i+1))){if(expected.length()>0){expected.tab();} expected.append(value, 9);}
					expected.tab().append(target);
					equal(file.nextLine(), expected.toString());
					expectedRows.append(FOLDS[i]).tab().append(row+1).append("\tg000").append(i+1)
							.tab().append(row==3 ? "control" : "candidate").tab().append(row==3 ? 1 : row==2 ? 4 : row+1)
							.tab().append(target).nl();
				}
				if(file.nextLine()!=null){throw new AssertionError("Unknown or extra trainer rows emitted");}
			}finally{if(file.close()){throw new AssertionError("Trainer table read failure");}}
		}
		final ByteFile rows=ByteFile.makeByteFile(root+"/rows.tsv", true);
		final ByteBuilder observed=new ByteBuilder();
		try{for(byte[] line=rows.nextLine(); line!=null; line=rows.nextLine()){observed.append(line).nl();}}
		finally{if(rows.close()){throw new AssertionError("Provenance read failure");}}
		equal(observed.toBytes(), expectedRows.toString());
		final ByteFile census=ByteFile.makeByteFile(root+"/census.tsv", true);
		final LineParser1 lp=new LineParser1('\t');
		try{
			census.nextLine();
			for(int i=0; i<3; i++){
				final byte[] line=census.nextLine();
				if(line==null){throw new AssertionError("Missing census row");}
				lp.set(line);
				final long[] expected={2, 1, 1, 1, 0, 1, 10, 2, 2, 7001+i};
				if(lp.terms()!=16){throw new AssertionError("Wrong census width");}
				for(int j=0; j<expected.length; j++){
					if(lp.parseLong(5+j)!=expected[j]){throw new AssertionError("Wrong class/sampling census");}
				}
			}
			if(census.nextLine()!=null){throw new AssertionError("Extra census rows");}
		}finally{if(census.close()){throw new AssertionError("Census read failure");}}
		System.err.println("SAME_K_TRAINING_FIXTURE_PASS rows=12 positive=9 negative=3 excluded=6");
	}

	/** Writes one tiny fixture using native I/O and observes writer errors. */
	private static void write(final String path, final ByteBuilder text){
		final ByteStreamWriter writer=new ByteStreamWriter(path, false, false, false);
		writer.start(); writer.print(text);
		if(writer.poisonAndWait()){throw new IllegalStateException("Fixture write failed: "+path);}
	}

	/** Byte equality detects even an unintended reformat of the producer's decimal fields. */
	private static void equal(final byte[] actual, final String expected){
		if(!Arrays.equals(actual, expected.getBytes(StandardCharsets.UTF_8))){throw new AssertionError("Fixture bytes differ from expected: "+expected);}
	}

	private static final String[] FOLDS={"train", "validation", "test"};
}
