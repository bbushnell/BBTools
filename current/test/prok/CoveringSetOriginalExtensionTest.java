package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.CoveringSet;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.RrnaResourceIO;
import shared.Shared;

/** Original prevalence must outrank the residual shortlist, with immediate
 * per-gene coverage, fixed prefix order and no batch overfill.
 * @author Brian Bushnell, Diona
 */
public final class CoveringSetOriginalExtensionTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, CoveringSetOriginalExtensionTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				check(parser.parse(arg, pair[0], pair.length>1 ? pair[1] : null), "Use out=fresh_directory");
			}
			check(parser.out1!=null && !parser.overwrite && !parser.append, "Use a fresh fixture directory");
			run(new File(parser.out1));
			pp.outstream.println("COVERING_SET_ORIGINAL_EXTENSION_FIXTURES_PASS");
		}finally{
			Shared.closeStream(pp.outstream);
		}
	}

	private static void run(final File root){
		check(root.mkdir(), "Fixture directory must be fresh");
		write(root, "initial.fa", ">first\nAA\n>second\nCC\n");
		write(root, "pool.fa", ">target1\nAACGTT\n>target2\nAACGTT\n"
			+">covered1\nAANNCCNNGT\n>covered2\nAANNCCNNGT\n>covered3\nAANNCCNNGT\n");
		invoke(root, "pool.fa", "legacy", "greedy", 1, "steps=1");
		invoke(root, "pool.fa", "original", "original", 1, "steps=1");
		invoke(root, "pool.fa", "parallel", "original", 4, "steps=4");
		keys(root, "legacy", "AA", "CC", "AC");
		keys(root, "original", "AA", "CC", "GT");
		keys(root, "parallel", "AA", "CC", "GT");
		summary(root, "original", 5, 2, 1, 0);
		summary(root, "parallel", 5, 2, 1, 0);

		write(root, "zero.fa", ">zero\nTTNNGGNNAG\n>covered\nAANNCC\n");
		invoke(root, "zero.fa", "zero", "original", 1, "steps=4");
		keys(root, "zero", "AA", "CC", "AG", "GG");
		summary(root, "zero", 2, 1, 2, 0);

		write(root, "repeat.fa", ">short_vocabulary\nGGGGGG\n>covered\nAANNCC\n");
		invoke(root, "repeat.fa", "repeat", "original", 1, "steps=4");
		keys(root, "repeat", "AA", "CC", "GG");
		summary(root, "repeat", 2, 1, 1, 1);

		write(root, "occurrences.fa", ">target\nAANNGGNNGT\n>repeat_copy\nAANNCCNNGGGGGG\n"
			+">copy1\nAANNCCNNGT\n>copy2\nAANNCCNNGT\n>copy3\nAANNCCNNGT\n>copy4\nAANNCCNNGT\n");
		invoke(root, "occurrences.fa", "occurrences", "original", 1, "steps=1");
		keys(root, "occurrences", "AA", "CC", "GG");
		summary(root, "occurrences", 6, 1, 1, 0);

		boolean rejected=false;
		try{
			new CoveringSet(new String[]{"in="+path(root, "pool.fa"), "out="+path(root, "invalid.fa"),
				"seedset="+path(root, "initial.fa"), "seedextension=original", "k=2", "target=1", "mincount=2"});
		}catch(IllegalArgumentException expected){
			rejected=true;
		}
		check(rejected, "Original extension cannot silently retain a conflicting residual-count filter");
	}

	private static void invoke(final File root, final String input, final String name,
			final String mode, final int threads, final String steps){
		CoveringSet.main(new String[]{"in="+path(root, input), "out="+path(root, name+".fa"),
			"summary="+path(root, name+".summary.tsv"), "seedset="+path(root, "initial.fa"),
			"seedextension="+mode, "k=2", "kdesign=2", "coveragek=2", "minhits=2", "mincount=1",
			"mininitial=0", "target=1", "rcomp=f", "ow=f", "pigz=f", "unpigz=f", "t="+threads, steps});
	}

	private static void keys(final File root, final String name, final String... expected){
		final ArrayList<String> actual=new ArrayList<String>();
		RrnaResourceIO.stream(path(root, name+".fa"), read->{
			actual.add(new String(read.bases, StandardCharsets.US_ASCII));
		});
		check(actual.equals(Arrays.asList(expected)), "Literal selected keys and fixed prefix order differ: "+name+" "+actual);
	}

	private static void summary(final File root, final String name, final int members,
			final int initialRemaining, final int added, final int remaining){
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(path(root, name+".summary.tsv"))){
			input.reader.nextLine();
			p.set(input.reader.nextLine());
			check(p.terms()==11 && p.parseInt(0)==members && p.parseInt(1)==2
				&& p.parseInt(2)==2+added && p.parseInt(3)==added && p.parseInt(7)==initialRemaining
				&& p.parseInt(8)==remaining && input.reader.nextLine()==null,
				"Original extension must fill only distinct missing slots and report irreducible residuals: "+name);
		}
	}

	private static void write(final File root, final String name, final String text){
		try(Output output=new Output(path(root, name))){
			output.writer.print(text);
		}
	}

	private static String path(final File root, final String name){
		return new File(root, name).getPath();
	}

	private static void check(final boolean value, final String why){
		if(!value){throw new AssertionError(why);}
	}
}
