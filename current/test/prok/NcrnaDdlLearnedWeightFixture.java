package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import parse.LineParser1;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.RankDualNcrnaWindows;
import structures.ByteBuilder;

/** Literal finite fit through the actual rank/report path, not a mocked weight.
 * Beta's one SEEN locus has two windows: A20 and A20-G20. Its background model
 * has per-window scores0/1 in both halves, so mean=.5, sd=.5 and learnedw=1.
 * Alpha's SEEN locus is in the judging group and has A20-G7: if leaked into
 * the short7 fit, its unequal half-scores would change w from1 to1/3.
 * @author Brian, Raiden
 */
final class NcrnaDdlLearnedWeightFixture{

	static void run(final String root){
		final RankSavedNcrnaWindowsTest.Fixture f=new RankSavedNcrnaWindowsTest.Fixture(root, "learned");
		RankDualNcrnaWindowsTest.coherent(f);
		rewrite(f.dir+"/panels.tsv", (p, r, b)->{
			if(alpha(p, 0)){r[8]="2";}
			fields(p, r, b);
		});
		rewrite(f.dir+"/selection.tsv", (p, r, b)->{
			if(alpha(p, 0) && p.parseInt(1)==0){
				r[6]="2";
				r[10]="0";
			}
			fields(p, r, b);
		});
		rewrite(f.dir+"/alpha.parents.tsv.gz", (p, r, b)->{if(!p.termEquals("other_subtype", 3)){fields(p, r, b);}});
		rewrite(f.dir+"/alpha.loci.tsv.gz", (p, r, b)->{
			if(p.termEquals("p1", 2)){
				r[6]="0";
				r[7]="1";
				r[8]="1";
				r[9]="RANKED";
			}
			fields(p, r, b);
		});
		rewrite(f.dir+"/alpha.windows.tsv.gz", (p, r, b)->{
			if(p.parseInt(2)==2){
				r[9]="126";
				r[11]="126";
				r[14]="27";
			}
			fields(p, r, b);
		});
		rewrite(f.dir+"/alpha.ranks.tsv.gz", (p, r, b)->{
			if(p.termEquals("p0", 5) || p.termEquals("p1", 5)){
				r[10]="1";
				r[11]="1";
				r[12]="1";
				r[13]="1";
				r[14]="0";
				r[15]="RANKED";
				if(p.termEquals("p1", 5)){r[7]="27";}
			}
			fields(p, r, b);
		});
		rewrite(f.dir+"/alpha.oracle.tsv.gz", (p, r, b)->{
			final int event=p.parseInt(2);
			if(event<3 || event>=5){
				r[7]=p.parseInt(5)==0 ? "12" : "0";
				r[8]=p.parseInt(5)==0 ? "1" : "0";
				if(event==2 && p.parseInt(5)==0){r[14]="1";}
			}
			fields(p, r, b);
		});
		rewrite(f.dir+"/alpha.links.tsv.gz", (p, r, b)->{
			if(p.termEquals("other_subtype", 6)){return;}
			if(p.parseInt(2)==2){
				r[8]="1";
				r[9]="true";
			}
			fields(p, r, b);
		});
		rewrite(f.dir+"/beta.inputs.tsv", (p, r, b)->{
			r[5]="2";
			fields(p, r, b);
		});
		rewrite(f.dir+"/beta.loci.tsv.gz", (p, r, b)->{
			r[4]="2";
			r[5]="2";
			r[7]="2";
			fields(p, r, b);
		});
		rewrite(f.dir+"/beta.windows.tsv.gz", (p, r, b)->{
			fields(p, r, b);
			r[2]="1";
			r[9]="39";
			r[11]="39";
			r[14]="40";
			fields(p, r, b);
		});
		rewrite(f.dir+"/beta.ranks.tsv.gz", (p, r, b)->{
			fields(p, r, b);
			r[2]="1";
			r[7]="40";
			r[13]="2";
			fields(p, r, b);
		});
		rewrite(f.dir+"/beta.oracle.tsv.gz", (p, r, b)->{
			fields(p, r, b);
			r[2]="1";
			r[7]="12";
			r[8]=p.parseInt(5)==0 ? "1" : "2";
			fields(p, r, b);
		});
		rewrite(f.dir+"/beta.links.tsv.gz", (p, r, b)->{
			fields(p, r, b);
			r[2]="1";
			fields(p, r, b);
		});
		sequences(f.dir);
		RankDualNcrnaWindowsTest.baseline(f, "baseline.tsv.gz", false, true);
		RankDualNcrnaWindows.run(RankDualNcrnaWindowsTest.config(f, "finite", 2));
		verify(f.dir+"/finite");
	}

	private static void verify(final String prefix){
		final LineParser1 p=new LineParser1('\t');
		boolean weight=false, judge=false, excluded=false;
		try(Input in=new Input(prefix+".weights.tsv")){
			check(in.reader.nextLine()!=null, "Learned weights require a header");
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				if(p.parseInt(0)==1 && p.termEquals("FIT_A_JUDGE_B", 1)){
					check(p.parseInt(2)==1 && number(p, 3)==1 && number(p, 4)==.5 && number(p, 5)==.5
						&& number(p, 9)==1 && number(p, 10)==1 && number(p, 11)==1 && p.termEquals("MEASURED", 13),
						"The outside-group SEEN locus alone fits w=1; including the judging-group SEEN locus would give short7 w=1/3");
					weight=true;
				}
			}
		}
		try(Input in=new Input(prefix+".summary.tsv")){
			check(in.reader.nextLine()!=null, "Matched held-out comparisons require a header");
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				if(!p.termEquals("g1_FIT_A_JUDGE_B", 0) || !p.termEquals("Q1_LOW", 2)){continue;}
				if(p.termEquals("NOT_SEEN", 1) && p.termEquals("JUDGING", 3)){
					check(p.parseInt(5)==1 && p.parseInt(6)==1 && p.parseInt(7)==1 && p.parseInt(11)==1
						&& p.parseInt(15)==1 && p.parseInt(19)==1, "Finite learned weight must reach the judging rank and all three matched reference sums");
					judge=true;
				}
				if(p.termEquals("SEEN", 1) && p.termEquals("EXCLUDED_SEEN_JUDGE_GROUP", 3)){
					check(p.parseInt(5)==1 && p.parseInt(11)==1, "Excluded SEEN reference sums must use that one-locus cohort, not the two-locus full SEEN population");
					excluded=true;
				}
			}
		}
		check(weight && judge && excluded, "The finite fit, judging rank and excluded-SEEN comparison must all be present");
	}

	private static void sequences(final String dir){
		final String alpha=dir+"/alpha.windows.fa.gz";
		final String original=preserve(alpha);
		try(Input in=new Input(original); Output out=new Output(alpha)){
			boolean second=false;
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				if(row.length>0 && row[0]=='>'){second=Arrays.equals(row, ">0_2".getBytes(StandardCharsets.US_ASCII));}
				out.writer.println(second && row[0]!='>' ? (A+"GGGGGGG").getBytes(StandardCharsets.US_ASCII) : row);
			}
		}
		final String beta=dir+"/beta.windows.fa.gz";
		preserve(beta);
		try(Output out=new Output(beta)){out.writer.print((">0_0\n"+A+"\n>0_1\n"+A+G+"\n").getBytes(StandardCharsets.US_ASCII));}
	}

	private interface Edit{void apply(LineParser1 p, String[] replacement, ByteBuilder out);}
	private static void rewrite(final String path, final Edit edit){
		final String original=preserve(path);
		final LineParser1 p=new LineParser1('\t');
		final String[] replacement=new String[16];
		final ByteBuilder b=new ByteBuilder();
		try(Input in=new Input(original); Output out=new Output(path)){
			out.writer.println(in.reader.nextLine());
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				Arrays.fill(replacement, null);
				b.clear();
				edit.apply(p, replacement, b);
				out.writer.print(b);
			}
		}
	}
	private static String preserve(final String path){
		final String original=path.endsWith(".gz") ? path.substring(0, path.length()-3)+".learned_original.gz" : path+".learned_original";
		check(!new File(original).exists() && new File(path).renameTo(new File(original)), "Preserve each original synthetic table before changing its literal evidence");
		return original;
	}
	private static void fields(final LineParser1 p, final String[] replacement, final ByteBuilder b){
		for(int field=0; field<p.terms(); field++){
			if(field>0){b.tab();}
			if(replacement[field]!=null){b.append(replacement[field]);}
			else{
				p.setBounds(field);
				b.append(p.line(), p.a(), p.b()-p.a());
			}
		}
		b.nl();
	}
	private static boolean alpha(final LineParser1 p, final int field){return Arrays.equals(p.parseByteArray(field), ALPHA);}
	/** A few fixture assertions read scientific notation; the native fast parser
	 * intentionally does not support exponents (Parse.parseDouble:499). */
	private static double number(final LineParser1 p, final int field){return Double.parseDouble(p.parseString(field));}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
	private static final String A="AAAAAAAAAAAAAAAAAAAA", G="GGGGGGGGGGGGGGGGGGGG";
	private static final byte[] ALPHA="α-panel".getBytes(StandardCharsets.UTF_8);
}
