package test.prok;

import java.io.File;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.Locale;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.TrnaConsensusBuilder;
import shared.Shared;

/** Unknown source taxonomy must preserve the assignment and every scientific
 * output byte. Two equal originals remain two members of their one cluster.
 * @author Brian Bushnell, Raiden
 */
public final class TrnaConsensusAssignmentsTaxidTest{

	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, TrnaConsensusAssignmentsTaxidTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown fixture argument: "+arg);
			}
			check(parser.out1!=null && !parser.overwrite && !parser.append && parser.in1==null,
				"Supply out=fresh_fixture_directory");
			final File root=new File(parser.out1).getAbsoluteFile();
			check(root.mkdir(), "The fixture directory must be fresh: "+root);
			final File input=new File(root, "members.fa");
			try(Output out=new Output(input.getPath())){
				out.writer.print(">known tid_123\n"+DNA+"\n>unknown\n"+DNA+"\n");
			}
			build(input, root, "without", false);
			build(input, root, "with", true);
			for(final String suffix:new String[]{".fa", ".hbm"}){
				final byte[] without=Files.readAllBytes(new File(root, "without"+suffix).toPath());
				final byte[] with=Files.readAllBytes(new File(root, "with"+suffix).toPath());
				check(without.length>0 && Arrays.equals(without, with),
					"Writing unknown-taxID assignments must not change consensus/HBM bytes: "+suffix);
			}
			final LineParser1 row=new LineParser1('\t');
			int count=0;
			boolean known=false, unknown=false;
			try(Input in=new Input(new File(root, "assignments.tsv").getPath())){
				final byte[] header=in.reader.nextLine();
				check(header!=null && header[0]=='#', "Assignment header is required");
				for(byte[] line; (line=in.reader.nextLine())!=null;){
					row.set(line);
					check(row.terms()==6 && row.termEquals("kept", 5), "Every original needs its intact six-column kept assignment");
					if(row.termEquals("known tid_123", 0)){
						check(!known && row.parseInt(1)==123, "Known taxonomy must be unchanged");
						known=true;
					}else if(row.termEquals("unknown", 0)){
						check(!unknown && row.parseInt(1)==-1, "Unknown taxonomy must remain -1 and retain its original identity");
						unknown=true;
					}else{throw new AssertionError("Unexpected assignment identity");}
					count++;
				}
			}
			check(count==2 && known && unknown, "Both original votes require exactly one assignment");
			pp.outstream.println("TRNA_CONSENSUS_ASSIGNMENTS_TAXID_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Uses the unchanged public builder path and same DNA/configuration in both runs. */
	private static void build(final File input, final File root, final String name, final boolean assignments){
		final java.util.ArrayList<String> args=new java.util.ArrayList<String>();
		args.add("in="+input.getPath());
		args.add("out="+new File(root, name+".fa"));
		args.add("outmodel="+new File(root, name+".hbm"));
		args.add("family=cobalamin");
		args.add("cluster=f");
		args.add("mincluster=1");
		args.add("mingroup=1");
		args.add("passes=2");
		args.add("endtrim=0");
		args.add("trd=f");
		args.add("t=1");
		args.add("overwrite=f");
		if(assignments){args.add("outassignments="+new File(root, "assignments.tsv"));}
		TrnaConsensusBuilder.main(args.toArray(new String[0]));
	}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
	private static final String DNA="ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
}
