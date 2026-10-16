package assemble;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import structures.ByteBuilder;

/** Routing/schema fixtures, not evidence of biological train/validation separation. @author Fischl */
public final class LiveBranchTrainingCorpusTest {

	/** Generates an isolated fixture or verifies exact partitioned rows against its native inputs. */
	public static void main(final String[] args){
		if(args.length==3 && args[0].equals("generate")){generate(args[1], args[2]);}
		else if(args.length==2 && args[0].equals("check")){checkOutput(args[1]);}
		else{throw new IllegalArgumentException("generate NEW_directory variant, or check fixture_directory");}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Uses the actual label exporter to exercise both sampled TRAIN and unsampled VALID contracts. */
	private static void generate(final String path, final String variant){
		final String root=new File(path).getAbsolutePath();
		if(!new File(root).mkdir()){throw new IllegalArgumentException("Fixture directory must be new.");}
		final ByteBuilder manifest=new ByteBuilder(LiveBranchTrainingCorpus.MANIFEST_HEADER).nl();
		for(int i=0; i<2; i++){
			final String fold=i==0 ? "train" : "validation", run=root+"/"+fold;
			LiveBranchCorpusTest.main(new String[]{"generate_sampling", run, "base"});
			final double fraction=i==0 || variant.equals("sampled_validation") ? 0.1 : 1;
			LiveBranchCorpus.main(new String[]{"in="+run+"/trace.tsv", "ref="+run+"/reference.fa", "out="+run+"/labeled",
					"controlfraction="+fraction, "controlseed=123", "t=1"});
			write(run+"/exit_code.txt", "0\n");
			manifest.append(i==0 ? "g0001" : "g0002").tab().append(fold).append("\t24\t30\t123\tnew\t")
					.append(run).tab().append(run+"/labeled").tab().append(run+"/unused.fq").tab().append(run+"/reference.fa").nl();
		}
		write(root+"/selected.tsv", manifest.toString());
		final String train=root+"/train/labeled/";
		if(variant.equals("duplicate_manifest")){appendLine(root+"/selected.tsv", 1);}
		else if(variant.equals("test_fold")){replace(root+"/selected.tsv", 2, "g0002\tvalidation\t", "g0002\ttest\t");}
		else if(variant.equals("wrong_k")){replace(root+"/selected.tsv", 1, "\t24\t30\t", "\t25\t30\t");}
		else if(variant.equals("wrong_seed")){replace(root+"/selected.tsv", 1, "\t123\t", "\t124\t");}
		else if(variant.equals("missing_genome")){removeLast(root+"/selected.tsv");}
		else if(variant.equals("wrong_header")){replace(train+"vectors.tsv", 0, "#dims\t12\t1", "#dims\t13\t1");}
		else if(variant.equals("bad_target")){
			final ArrayList<String> rows=lines(train+"rows.tsv");
			final String row=rows.get(1);
			rows.set(1, row.substring(0, row.length()-1)+(row.endsWith("1") ? "0" : "1"));
			writeLines(train+"rows.tsv", rows);
		}else if(variant.equals("duplicate_word")){
			final ArrayList<String> rows=lines(train+"rows.tsv");
			final String[] a=rows.get(1).split("\t"), b=rows.get(2).split("\t");
			rows.set(2, b[0]+"\t"+a[1]+"\t"+b[2]+"\t"+b[3]); writeLines(train+"rows.tsv", rows);
		}else if(variant.equals("extra_vector")){appendLine(train+"vectors.tsv", 2);}
		else if(variant.equals("truncated_rows")){removeLast(train+"rows.tsv");}
		else if(variant.equals("nonfinite") || variant.equals("wrong_width")){
			final ArrayList<String> rows=lines(train+"vectors.tsv");
			final String row=rows.get(2);
			rows.set(2, (variant.equals("nonfinite") ? "NaN" : "")+row.substring(row.indexOf('\t')+(variant.equals("wrong_width") ? 1 : 0)));
			writeLines(train+"vectors.tsv", rows);
		}else if(variant.equals("bad_completion")){replace(train+"complete.txt", 0, "controls=1 ", "controls=2 ");}
		else if(variant.equals("bad_census")){replace(train+"census.tsv", 1, "grey\t9\t0\t0", "grey\t8\t0\t0");}
		else if(variant.equals("unexpected_role")){replace(train+"rows.tsv", 1, "\tgrey\t", "\thard\t");}
	}

	/** Checks every emitted value line against the correct native input population. */
	private static void checkOutput(final String root){
		checkPartition(root, "train", "train", null, 10);
		checkPartition(root, "validation", "validation", "grey", 9);
		checkPartition(root, "validation_controls", "validation", "control", 20);
		final ByteFile rows=ByteFile.makeByteFile(root+"/result/rows.tsv", true);
		final LineParser1 fields=new LineParser1('\t');
		final long[] counts=new long[3];
		try{
			check("partition\trow\tgenome\tk\tword\trole\ttarget".equals(text(rows.nextLine())), "Wrong merged provenance header.");
			for(byte[] line=rows.nextLine(); line!=null; line=rows.nextLine()){
				fields.set(line);
				final int partition=fields.termEquals("train", 0) ? 0 : fields.termEquals("validation", 0) ? 1 :
						fields.termEquals("validation_controls", 0) ? 2 : -1;
				check(partition>=0 && fields.terms()==7 && fields.parseLong(1)==++counts[partition], "Incorrect partition ordinal.");
				check(fields.termEquals(partition==0 ? "g0001" : "g0002", 2) && fields.parseInt(3)==24, "Genome leaked into another fold.");
				if(partition>0){check(fields.termEquals(partition==1 ? "grey" : "control", 5), "Validation checkpoint population was mixed.");}
			}
		}finally{if(rows.close()){throw new IllegalStateException("Cannot read merged provenance.");}}
		check(counts[0]==10 && counts[1]==9 && counts[2]==20, "Wrong merged population counts.");
		System.err.println("LIVE_BRANCH_TRAINING_TEST_PASS train=10 validation=9 validation_controls=20");
	}

	/** Native vectors must remain byte-identical; only routing changes. */
	private static void checkPartition(final String root, final String partition, final String source, final String role, final int expected){
		final ByteFile values=ByteFile.makeByteFile(root+"/"+source+"/labeled/vectors.tsv", true);
		final ByteFile rows=ByteFile.makeByteFile(root+"/"+source+"/labeled/rows.tsv", true);
		final ByteFile output=ByteFile.makeByteFile(root+"/result/"+partition+".tsv", true);
		final LineParser1 fields=new LineParser1('\t');
		int count=0;
		try{
			check(values.nextLine()!=null && values.nextLine()!=null && rows.nextLine()!=null, "Missing source headers.");
			check("#dims\t12\t1".equals(text(output.nextLine())) && LiveBranchTrainingCorpus.SCHEMA.equals(text(output.nextLine())), "Wrong merged dimensions/schema.");
			for(byte[] row=rows.nextLine(); row!=null; row=rows.nextLine()){
				fields.set(row);
				final byte[] vector=values.nextLine();
				check(vector!=null, "Missing native vector.");
				if(role==null || fields.termEquals(role, 2)){
					check(Arrays.equals(vector, output.nextLine()), "Partition values differ from their native source."); count++;
				}
			}
			check(values.nextLine()==null && output.nextLine()==null && count==expected, "Partition omitted or added vectors.");
		}finally{if(values.close() | rows.close() | output.close()){throw new IllegalStateException("Cannot compare partition inputs.");}}
	}

	/** Small fixture records only; production conversion does not retain row Strings. */
	private static ArrayList<String> lines(final String path){
		final ArrayList<String> lines=new ArrayList<String>();
		final ByteFile file=ByteFile.makeByteFile(path, true);
		try{for(byte[] row=file.nextLine(); row!=null; row=file.nextLine()){lines.add(text(row));}}
		finally{if(file.close()){throw new IllegalStateException("Cannot read fixture.");}}
		return lines;
	}

	/** Mutation is confined to a newly generated negative-test directory. */
	private static void replace(final String path, final int row, final String before, final String after){
		final ArrayList<String> lines=lines(path);
		check(lines.get(row).contains(before), "Negative fixture mutation did not match.");
		lines.set(row, lines.get(row).replace(before, after)); writeLines(path, lines);
	}

	/** Deliberately duplicates one generated input row. */
	private static void appendLine(final String path, final int index){
		final ArrayList<String> lines=lines(path); lines.add(lines.get(index)); writeLines(path, lines);
	}

	/** Deliberately truncates a generated input. */
	private static void removeLast(final String path){
		final ArrayList<String> lines=lines(path); lines.remove(lines.size()-1); writeLines(path, lines);
	}

	/** Keeps fixture rewriting on the native byte writer. */
	private static void writeLines(final String path, final ArrayList<String> lines){
		final ByteBuilder out=new ByteBuilder();
		for(String line : lines){out.append(line).nl();}
		write(path, out.toString());
	}

	/** Overwriting is permitted only for deliberately malformed scratch fixtures. */
	private static void write(final String path, final String text){
		final ByteStreamWriter writer=new ByteStreamWriter(path, true, false, false);
		writer.start(); writer.print(text);
		if(writer.poisonAndWait()){throw new IllegalStateException("Cannot write test fixture.");}
	}

	private static String text(final byte[] line){return line==null ? null : new String(line, StandardCharsets.US_ASCII);}
	private static void check(final boolean ok, final String message){if(!ok){throw new AssertionError(message);}}
}
