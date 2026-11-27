package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Locale;

import map.LongIntMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.ReportNcrnaSeedPositions;
import prok.ReportNcrnaSeedPositions.Key;
import prok.ReportNcrnaSeedPositions.Scratch;
import shared.Shared;

/** Literal reference-only oracles; no external genes or caller execution.
 * @author Brian, Diona, Raiden
 */
public final class ReportNcrnaSeedPositionsTest{
	/** Creates a fresh fixture directory and owns all literal, external-input-free checks. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, ReportNcrnaSeedPositionsTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String a=pair[0].toLowerCase(Locale.ROOT), b=pair.length>1?pair[1]:null;
				check(parser.parse(arg, a, b), "Use native out=/thread/JVM options");
			}
			check(parser.out1!=null && parser.in1==null && !parser.overwrite && !parser.append, "Fixtures need a fresh directory");
			final String root=new File(parser.out1).getAbsolutePath();
			check(new File(root).mkdir(), "Preserve earlier fixture outputs");
			quantiles();
			positions();
			pipeline(root);
			pp.outstream.println("NCRNA_SEED_POSITION_FIXTURES_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Checks the small-sample p99.9 boundary and the declared nearest-rank convention. */
	private static void quantiles(){
		final int[] a=new int[999], b=new int[1000];
		a[998]=100;
		b[999]=100;
		check(ReportNcrnaSeedPositions.percentile(a, 999, 999, 1000)==100
			&& ReportNcrnaSeedPositions.percentile(b, 1000, 999, 1000)==0,
			"Nearest-rankp999 uses the maximum below1000 samples; do not interpolate or silently usep99.95");
		check(ReportNcrnaSeedPositions.percentile(new int[]{17, 35, 41, 41}, 4, 1, 2)==35, "Median is the stated lower nearest rank for even samples");
		boolean rejected=false;
		try{
			ReportNcrnaSeedPositions.rank(0, 999, 1000);
		}catch(IllegalArgumentException expected){
			rejected=true;
		}
		check(rejected, "Empty-key caps must stay unavailable, never an invented zero percentile");
	}
	/** Literal zero-based oracles separate seed-edge caps from native center-based blanket padding. */
	private static void positions(){
		final LongIntMap index=new LongIntMap(8);
		index.put(encode(A), 0);
		index.put(encode(B), 1);
		index.put(encode(C), 2);
		final Scratch scratch=new Scratch(3);
		final String sequence="GG"+A+"CCC"+B+"AA";
		scratch.scan(bytes(sequence.toLowerCase(Locale.ROOT)), index);
		check(scratch.hitKeys.size==2 && scratch.distinct.size==2 && scratch.starts.get(0)==2 && scratch.starts.get(1)==22,
			"Literal two17-mers begin at2and22; lowercase uses the native forward lookup");
		scratch.scan(bytes(A.substring(0, 8)+"N"+A.substring(8)), index);
		check(scratch.hitKeys.size==0, "A seed cannot bridge an ambiguous base");
		scratch.scan(bytes("AAAAAAAAAAAAAAAAA"), index);
		check(scratch.hitKeys.size==0, "The reverse complement of the selected T17 key is not a forward hit");
		scratch.scan(bytes(A+"N"+A), index);
		check(scratch.hitKeys.size==2 && scratch.distinct.size==1 && scratch.starts.get(0)==0 && scratch.starts.get(1)==18,
			"Repeated-key occurrences remain separate observations but one distinct key");
		final Key[] keys={new Key(), new Key(), new Key()};
		keys[0].leftCap=18;
		keys[0].rightCap=22;
		keys[1].leftCap=22;
		keys[1].rightCap=2;
		scratch.proposals(keys, 0, 24);
		check(scratch.combineStart(0)==-18 && scratch.combineStop(0)==56 && scratch.combineStart(1)==-18 && scratch.combineStop(1)==56,
			"Repeated occurrences first form one key envelope; they are not two independent intersection votes");
		scratch.scan(bytes(sequence), index);
		scratch.proposals(keys, 0, 24);
		check(scratch.combineStart(0)==-16 && scratch.combineStop(0)==40 && scratch.combineStart(1)==0 && scratch.combineStop(1)==40,
			"Per-key union and intersection retain exact seed-edge geometry without genomic clipping");
		scratch.proposals(keys, 1, 24);
		check(scratch.combineStart(0)==-14 && scratch.combineStop(0)==71 && scratch.combineStart(1)==6 && scratch.combineStop(1)==51,
			"Native blanket geometry uses centers; its intersection can miss the true5prime end");
		scratch.proposals(keys, 2, 24);
		check(scratch.hitLengths.get(0)==918 && scratch.containedHits==2 && scratch.combineStart(0)==-440 && scratch.combineStop(0)==497,
			"The stated legacy comparison is918bases per hit, including native inclusive endpoints");
		keys[0].leftCap=0;
		keys[0].rightCap=0;
		keys[1].leftCap=0;
		keys[1].rightCap=0;
		scratch.proposals(keys, 0, 24);
		check(scratch.combineStop(1)<scratch.combineStart(1) && scratch.containedHits==0, "Disjoint key intervals must remain an explicit empty intersection and failed per-hit containment");
	}
	/** Runs four sense-oriented genes; minus annotations do not reverse their supplied bases. */
	private static void pipeline(final String root){
		final String genes=root+"/genes.fa", seeds=root+"/seeds.fa", prefix=root+"/result";
		try(Output out=new Output(genes)){
			out.writer.println(">plus strand=+\nGG"+A+"CCC"+B+"AA\n>minus strand=-\nGG"+A+"CCC"+B+"AA\n>repeat strand=+\n"+A+"N"+A+"\n>zero strand=-\nNNNNNNNNNNNNNNNNN");
		}
		try(Output out=new Output(seeds)){out.writer.println(">a\n"+A+"\n>b\n"+B+"\n>unseen\n"+C);}
		ReportNcrnaSeedPositions.main(new String[]{"in="+genes, "seeds="+seeds, "member=fixture", "expectedgenes=4", "expectedkeys=3", "orientation=sense", "out="+prefix, "t=1"});
		final LineParser1 p=new LineParser1('\t');
		int rows=0;
		try(Input in=new Input(prefix+".keys.tsv")){
			in.reader.nextLine();
			for(byte[] line; (line=in.reader.nextLine())!=null; ){
				p.set(line);
				rows++;
				if(p.termEquals(A, 1)){
					check(p.parseInt(2)==4 && p.parseInt(3)==3 && p.parseInt(4)==18 && p.parseInt(5)==22 && p.parseInt(6)==57,
						"Repeated A occurrences/counts/caps must retain literal reference weights");
				}else if(p.termEquals(B, 1)){
					check(p.parseInt(2)==2 && p.parseInt(3)==2 && p.parseInt(4)==22 && p.parseInt(5)==2 && p.parseInt(6)==41,
						"Minus headers do not reverse already sense-oriented genes");
				}else{
					check(p.termEquals(C, 1) && p.parseInt(2)==0 && p.parseInt(3)==0 && p.termEquals("NA", 4) && p.termEquals("NO_OCCURRENCES", 7),
						"Unobserved seed axis rows must remain explicit");
				}
			}
		}
		check(rows==3, "Every selected key has one result row");
		try(Input in=new Input(prefix+".lengths.tsv")){
			in.reader.nextLine();
			p.set(in.reader.nextLine());
			check(p.parseInt(1)==4 && p.parseInt(2)==17 && p.parseInt(3)==35 && p.parseInt(5)==41 && p.parseInt(7)==24
				&& p.parseInt(13)==1 && p.parseInt(14)==1 && p.parseInt(15)==2, "Gene population,nearest-rank lengths,padding and distinct-key strata remain separate");
		}
		int references=0, absent=0, plus=0, minus=0;
		try(Input in=new Input(prefix+".references.tsv")){
			in.reader.nextLine();
			for(byte[] line; (line=in.reader.nextLine())!=null; ){
				p.set(line);
				references++;
				if(p.termEquals("zero", 1)){
					absent++;
					check(p.termEquals("false", 7) && p.termEquals("NA", 8), "Zero-key genes have absent,not zero-sized,windows");
				}
				if(p.termEquals("plus", 1)){plus++;}
				if(p.termEquals("minus", 1)){minus++;}
			}
		}
		check(references==24 && absent==6 && plus==6 && minus==6, "All four references retain all six method/rule rows");
		int hitRows=0, allRows=0, zeroRows=0;
		try(Input in=new Input(prefix+".windows.tsv")){
			in.reader.nextLine();
			for(byte[] line; (line=in.reader.nextLine())!=null; ){
				p.set(line);
				if(p.termEquals("PER_HIT", 2)){
					hitRows++;
					check(p.parseInt(4)==6 && p.parseInt(7)==(p.termEquals("MEMBER_LENGTH_P999_CAP", 1)?3:6)
						&& p.parseInt(8)==(p.termEquals("MEMBER_LENGTH_P999_CAP", 1)?3:0), "Per-hit containment preserves its six occurrence denominator separately from reference windows");
					check(p.parseDouble(14)==(p.termEquals("MEMBER_LENGTH_P999_CAP", 1)?0.5:0.0), "Per-hit failure fractions use all six individual proposals");
				}else if(p.termEquals("ALL", 3)){
					allRows++;
					final boolean shortened=p.termEquals("MEMBER_LENGTH_P999_CAP", 1) && p.termEquals("INTERSECTION", 2);
					check(p.parseInt(4)==4 && p.parseInt(5)==1 && p.parseInt(7)==(shortened?1:3)
						&& p.parseDouble(15)==(shortened?0.75:0.25), "All-reference failure fractions keep the one no-hit reference in their denominator");
				}else if(p.termEquals("ZERO_KEYS", 3)){
					zeroRows++;
					check(p.termEquals("NA", 14) && p.parseDouble(15)==1.0, "Absent windows have unavailable conditional containment, but all-reference failure is one");
				}
			}
		}
		check(hitRows==3 && allRows==6 && zeroRows==6, "Each method has per-hit and both combined-reference summaries with explicit zero strata");
		lengthsOnly(root, genes);
	}
	/** Independent length-only output must retain all records without seed or window outputs. */
	private static void lengthsOnly(final String root, final String genes){
		final String prefix=root+"/lengths_only";
		ReportNcrnaSeedPositions.main(new String[]{"mode=lengths", "in="+genes, "member=fixture", "expectedgenes=4", "out="+prefix, "t=1"});
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(prefix+".lengths.tsv")){
			in.reader.nextLine();
			p.set(in.reader.nextLine());
			check(p.terms()==10 && p.parseInt(1)==4 && p.parseInt(2)==17 && p.parseInt(3)==35
				&& p.parseInt(4)==41 && p.parseInt(5)==41 && p.parseInt(6)==41
				&& p.parseInt(7)==24 && p.parseInt(8)==24 && in.reader.nextLine()==null,
				"One-pass reference lengths preserve four literal genes, exact nearest ranks and both minus17 caps");
		}
		check(!new File(prefix+".keys.tsv").exists() && !new File(prefix+".windows.tsv").exists()
			&& !new File(prefix+".references.tsv").exists() && !new File(prefix+".window_histogram.tsv").exists(),
			"Length-only mode must not emit unmeasured seed/window artifacts");
		final String duplicateInput=root+"/duplicate_names.fa", observed=root+"/observed_lengths";
		try(Output out=new Output(duplicateInput)){
			out.writer.println(">same\n"+A+"\n>same\n"+A+"N"+A);
		}
		ReportNcrnaSeedPositions.main(new String[]{"mode=lengths", "in="+duplicateInput, "member=duplicates", "out="+observed, "t=1"});
		try(Input in=new Input(observed+".lengths.tsv")){
			in.reader.nextLine();
			p.set(in.reader.nextLine());
			check(p.parseInt(1)==2 && p.parseInt(2)==17 && p.parseInt(3)==17 && p.parseInt(5)==35
				&& p.parseInt(6)==35 && p.parseInt(7)==18 && p.parseInt(8)==18,
				"Optional expected count and repeated names preserve two physical records of lengths17and35");
		}
	}
	/** Encodes only the three literal 17-base fixture keys, independently of the reporter scanner. */
	private static long encode(final String word){
		check(word.length()==17, "Literal fixture seeds are exactly17bases");
		long result=0;
		for(int i=0; i<word.length(); i++){
			final int digit="ACGT".indexOf(word.charAt(i));
			check(digit>=0, "Literal DNA alphabet");
			result=(result<<2)|digit;
		}
		return result;
	}
	/** Returns ASCII bases for a literal fixture; not an input-data parser. */
	private static byte[] bytes(final String text){return text.getBytes(StandardCharsets.US_ASCII);}
	/** Fails the fixture loudly even when assertions were accidentally disabled. */
	private static void check(final boolean ok, final String why){if(!ok){throw new AssertionError(why);}}

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	private static final String A="ACGTACGTACGTACGTA", B="TTTTTTTTTTTTTTTTT", C="CCCCCCCCCCCCCCCCC";
}
