package test.prok;

import java.io.File;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.ReportNcrnaSeedSupport;
import shared.Shared;

/** Literal occurrences differ from gene-copy support; a batch-only extra key
 * earns no first-two coverage credit. Repeated record IDs remain separate copies.
 * @author Brian Bushnell, Diona
 */
public final class ReportNcrnaSeedSupportTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, ReportNcrnaSeedSupportTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				check(parser.parse(arg, pair[0], pair.length>1 ? pair[1] : null), "Use out=fresh_directory");
			}
			check(parser.out1!=null && !parser.overwrite && !parser.append, "Use a fresh fixture directory");
			run(new File(parser.out1));
			pp.outstream.println("NCRNA_SEED_SUPPORT_FIXTURES_PASS");
		}finally{
			Shared.closeStream(pp.outstream);
		}
	}

	private static void run(final File root){
		check(root.mkdir(), "Fixture directory must be fresh");
		seeds(root, "initial.fa", A, C);
		seeds(root, "after.fa", A, C, G, T, X);
		try(Output output=new Output(path(root, "genes.fa"))){
			output.writer.print(">same\n"+A+'N'+G+"\n>one\n"+C+'N'+T+"\n>same\n"+A+'N'+C+'N'+G+'N'+G
				+"\n>zero\n"+G+'N'+T+"\n>already_covered\n"+A+'N'+C+'N'+X+'\n');
		}
		invoke(root, "after.fa", "result", 5, 5);
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(path(root, "result.summary.tsv"))){
			input.reader.nextLine();
			p.set(input.reader.nextLine());
			final long[] expected={5, 2, 5, 3, 2, 1, 1, 4, 4, 0, 1, 2, 2, 0, 0, 5};
			check(p.terms()==17, "Summary schema must keep support and missing-slot counts separate");
			for(int i=0; i<expected.length; i++){
				check(p.parseLong(i)==expected[i], "Wrong literal support/coverage count at column "+i);
			}
			check(input.reader.nextLine()==null, "Exactly one complete original-pool summary is required");
		}
		int rows=0;
		try(Input input=new Input(path(root, "result.keys.tsv"))){
			input.reader.nextLine();
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				rows++;
				check(p.parseInt(0)==rows, "Selection order must remain explicit");
				if(p.termEquals(G, 1)){
					check(p.parseLong(3)==4 && p.parseLong(4)==3 && p.parseLong(5)==2,
						"Four original occurrences in three copies fill two missing slots");
				}
				if(p.termEquals(X, 1)){
					check(p.termEquals("false", 2) && p.parseLong(3)==1 && p.parseLong(4)==1 && p.parseLong(5)==0,
						"A key added to an already-covered copy has no missing-slot credit");
				}
			}
		}
		check(rows==5, "All fixed and added keys must be reported");

		seeds(root, "reordered.fa", C, A, G, T, X);
		reject(()->invoke(root, "reordered.fa", "bad_order", 5, 5), "The fixed prefix cannot be reordered");
		reject(()->invoke(root, "after.fa", "bad_count", 5, 6), "Declared gene copies must match the complete input");
		try(Output output=new Output(path(root, "duplicate.fa"))){
			for(int i=0; i<50000; i++){output.writer.print(">duplicate\n"+G+'\n');}
		}
		reject(()->invoke(root, "duplicate.fa", "bad_duplicate", 50000, 5),
			"Duplicate ordered keys must fail after safely draining queued native input");
	}

	private static void invoke(final File root, final String seeds, final String out,
			final int expectedKeys, final int expectedGenes){
		ReportNcrnaSeedSupport.main(new String[]{"in="+path(root, seeds), "seedset="+path(root, "initial.fa"),
			"genes="+path(root, "genes.fa"), "counts="+expectedGenes, "expectedinitial=2",
			"expectedkeys="+expectedKeys, "out="+path(root, out), "t=1"});
	}

	private static void seeds(final File root, final String path, final String... keys){
		try(Output output=new Output(path(root, path))){
			for(int i=0; i<keys.length; i++){output.writer.print(">key"+i+'\n'+keys[i]+'\n');}
		}
	}

	private static String path(final File root, final String name){
		return new File(root, name).getPath();
	}

	private static void reject(final Runnable action, final String why){
		boolean failed=false;
		try{
			action.run();
		}catch(IllegalArgumentException expected){
			failed=true;
		}
		check(failed, why);
	}

	private static void check(final boolean value, final String why){
		if(!value){throw new AssertionError(why);}
	}

	private static final String A="AAAAAAAAAAAAAAAAA";
	private static final String C="CCCCCCCCCCCCCCCCC";
	private static final String G="GGGGGGGGGGGGGGGGG";
	private static final String T="TTTTTTTTTTTTTTTTT";
	private static final String X="ACACACACACACACACA";
}
