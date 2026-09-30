package prok;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;

/** Planted anchor, end geometry, shifted inference, sparse roundtrip and unequal
 * positional-denominator fixtures for the per-consensus endpoint score table.
 * @author Ganyu
 */
public final class RrnaPositionalKmerTableTest {
	public static void main(String[] args) throws Exception{
		final Path root=Files.createTempDirectory("rrna_positional_");final byte[] bases=new byte[180];Arrays.fill(bases,(byte)'A');
		final byte[] first="CGTACG".getBytes(StandardCharsets.US_ASCII),last="TGCATG".getBytes(StandardCharsets.US_ASCII);
		System.arraycopy(first,0,bases,47,6);System.arraycopy(last,0,bases,125,6);final String seq=new String(bases,StandardCharsets.US_ASCII);
		final StringBuilder fasta=new StringBuilder();for(int i=0;i<100;i++){fasta.append(">r").append(i).append(" model=test lflank=50 rflank=50\n").append(seq).append('\n');}
		final Path input=root.resolve("training.fa");Files.write(input,fasta.toString().getBytes(StandardCharsets.US_ASCII));final String prefix=root.resolve("table").toString();
		RrnaPositionalKmerTable.main(new String[]{"mode=build","in="+input,"model=test","expected=100","k=6","r=12","startoffset=-3","stopoffset=-4","pseudocount=0.5","outprefix="+prefix});
		final RrnaPositionalKmerTable.Table start=RrnaPositionalKmerTable.load(prefix+".5prime.tsv.gz","test"),stop=RrnaPositionalKmerTable.load(prefix+".3prime.tsv.gz","test");
		check(start.offset==-3&&stop.offset==-4&&start.records==100&&stop.records==100,"Both signed offsets and member counts must roundtrip without interpreting stop as an exclusive coordinate");
		final float[] values=new float[25];start.scoreWindow(bases,47,values);check(best(values)==12&&values[12]>0,"Planted 5prime anchor must uniquely outrank every shifted position");
		stop.scoreWindow(bases,125,values);check(best(values)==12&&values[12]>0,"Planted3prime k-mer uses last-base endpoint plus its declared start offset");
		start.scoreWindow(bases,50,values);check(best(values)==9,"Moving the inferred anchor3bases right must move the correct peak3positions left; do not apply offset twice");
		start.scoreWindow(bases,0,values);for(int i=0;i<12;i++){check(Float.isNaN(values[i]),"Clipped candidate scores must be explicit missing values, not neutral zeros");}
		final int word=RrnaPositionalKmerTable.encode(first,0,6);final double expected=Math.log((100+.5)/.5);
		check(Math.abs(start.scoreKey(word)-expected)<1e-12,"Background excludes the anchor; with no off-anchor copies the score is ln((100+.5)/.5)");
		check(Math.abs(start.scoreKey(RrnaPositionalKmerTable.encode("CCCCCC".getBytes(StandardCharsets.US_ASCII),0,6)))<1e-12,"Unseen key has a neutral prior when all positional sample sizes agree");
		final Path candidates=root.resolve("candidates.fa"),scores=root.resolve("scores.tsv");Files.write(candidates,(">candidate model=test anchor=47\n"+seq+'\n').getBytes(StandardCharsets.US_ASCII));
		RrnaPositionalKmerTable.main(new String[]{"mode=score","in="+candidates,"model=test","expected=1","table="+prefix+".5prime.tsv.gz","out="+scores});
		final String[] lines=new String(Files.readAllBytes(scores),StandardCharsets.US_ASCII).trim().split("\n"),fields=lines[1].split("\t");check(lines.length==2&&fields.length==30&&fields[4].equals("25"),"CLI emits exactly25 scores plus5metadata columns per candidate");
		check(Float.parseFloat(fields[17])>0,"CLI position0 agrees with direct API anchor score");
		boolean failed=false;try{RrnaPositionalKmerTable.load(prefix+".5prime.tsv.gz","wrong");}catch(IllegalArgumentException e){failed=true;}check(failed,"Wrong consensus table cannot load silently");
		final Path bad=root.resolve("bad_table.tsv");final String text=RrnaEndpointTestSupport.read(prefix+".5prime.tsv.gz");Files.write(bad,text.replace("test\t5prime\t6\t12\t","test\t5prime\t6\t13\t").getBytes(StandardCharsets.US_ASCII));
		failed=false;try{RrnaPositionalKmerTable.load(bad.toString(),"test");}catch(IllegalArgumentException e){failed=true;}check(failed,"Wrong radius must fail rather than alter inference width");
		failed=false;try{RrnaPositionalKmerTable.checkModel("row model=other","test");}catch(IllegalArgumentException e){failed=true;}check(failed,"Mixed consensus assignments must fail");
		testMissingDenominators();
		testExponentRoundtrip(root,text,start);
		System.out.println("PASS RrnaPositionalKmerTableTest planted_anchor both_ends shifted_window background_excludes_anchor per_position_denominators missing_scores cli_roundtrip wrong_model_radius exponent_roundtrip");
	}
	static void testExponentRoundtrip(Path root,String text,RrnaPositionalKmerTable.Table original)throws Exception{
		final StringBuilder output=new StringBuilder();
		for(String line:text.split("\n")){
			final String[] fields=line.split("\t",-1);
			if(fields.length==29 && fields[0].length()==6){
				fields[3]=String.format(java.util.Locale.ROOT,"%.17e",Double.parseDouble(fields[3]));
				line=String.join("\t",fields);
			}else if(fields.length==7 && fields[0].equals("test")){
				fields[5]="5E-1";line=String.join("\t",fields);
			}
			output.append(line).append('\n');
		}
		final Path file=root.resolve("exponents.tsv");Files.write(file,output.toString().getBytes(StandardCharsets.US_ASCII));
		final RrnaPositionalKmerTable.Table loaded=RrnaPositionalKmerTable.load(file.toString(),"test");
		check(loaded.alpha==original.alpha,"Scientific notation pseudocount must preserve table normalization");
		for(int word=0;word<original.space;word++){check(loaded.scoreKey(word)==original.scoreKey(word),"Exponent serialization must preserve every reconstructed score");}
		final parse.LineParser1 p=new parse.LineParser1('\t');p.set("-8.445034036997808E-4\t1e+2".getBytes(StandardCharsets.US_ASCII));
		check(RrnaPositionalKmerTable.savedDouble(p,0)==-8.445034036997808E-4 && RrnaPositionalKmerTable.savedDouble(p,1)==100,
			"The observed real-table small negative exponent and lowercase positive exponent must both load");
	}
	static void testMissingDenominators(){
		final byte[] clean=new byte[90];Arrays.fill(clean,(byte)'A');final byte[] word="CGTACG".getBytes(StandardCharsets.US_ASCII);System.arraycopy(word,0,clean,30,6);final byte[] n=clean.clone();n[41]='N';
		final RrnaPositionalKmerTable.Table t=new RrnaPositionalKmerTable.Table("test","5prime",6,0,.5);t.add(clean,30);t.add(n,30);t.finish();
		final double anchor=2.5/2050,background=(6*.5/2049+18*.5/2050)/24,expected=Math.log(anchor/background);
		check(Math.abs(t.scoreKey(RrnaPositionalKmerTable.encode(word,0,6))-expected)<1e-12,"Mean background frequency must normalize each offset: six positions have one valid record, eighteen have two");
		for(int site=0;site<25;site++){check(t.valid[site]+t.clipped[site]+t.ambiguous[site]==2,"Missing-base accounting must conserve every training observation");}
		final RrnaPositionalKmerTable.Table empty=new RrnaPositionalKmerTable.Table("test","5prime",6,0,.5);final byte[] ambiguous=new byte[80];Arrays.fill(ambiguous,(byte)'N');empty.add(ambiguous,30);boolean failed=false;
		try{empty.finish();}catch(IllegalArgumentException e){failed=true;}check(failed,"Pseudocount-only data cannot pretend to be trained evidence");
	}
	static int best(float[] values){int best=-1;for(int i=0;i<values.length;i++){if(Float.isFinite(values[i])&&(best<0||values[i]>values[best])){best=i;}}check(best>=0,"Fixture must have usable scores");for(int i=0;i<values.length;i++){if(i!=best){check(values[i]<values[best],"Planted fixture requires a unique maximum");}}return best;}
	static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
}
