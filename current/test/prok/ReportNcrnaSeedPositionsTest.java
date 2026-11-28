package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Locale;
import java.util.ArrayList;
import java.util.Arrays;

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
import prok.NcrnaOverlapWindowBuilder;
import prok.NcrnaOverlapWindowBuilder.Caps;
import prok.TrnaCaller;
import prok.GeneModel;
import prok.FamilyResourceResolver;
import prok.FamilyResourceResolver.Overrides;
import prok.NcrnaFamilyLoader;
import prok.CobalaminCallerResources;
import prok.ProkObject;
import prok.Orf;
import prok.NcrnaFamily;
import prok.NcrnaScavenger;
import prok.GeneCaller;
import prok.CallGenes;
import map.LongHashSet;
import shared.Shared;
import structures.ByteBuilder;

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
			parameterBounds(root);
			combiningBounds(root);
			windowGeometry();
			memberAttribution();
			resourceBinding(root);
			parallelTrnaCounts();
			ordinaryWorker();
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
		Arrays.fill(a, 100);
		Arrays.fill(b, 100);
		a[998]=10000;
		b[999]=10000;
		check(ReportNcrnaSeedPositions.chainCeiling(a, a.length)==9983
			&& ReportNcrnaSeedPositions.chainCeiling(b, b.length)==83,
			"At999originals the maximum survives; at1000the single giant outlier lies above nearest-rank p99.9");
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
		bounds(root, genes, seeds);
	}
	/** Producer/consumer round trip keeps metadata distinct from per-key margins. */
	private static void bounds(final String root, final String genes, final String seeds){
		final String prefix=root+"/bounds";
		ReportNcrnaSeedPositions.main(new String[]{"mode=bounds", "metadata=t", "in="+genes, "seeds="+seeds,
			"member=fixture", "expectedgenes=4", "expectedkeys=3", "orientation=sense", "out="+prefix, "t=1"});
		final long[] keys={encode(A), encode(B), encode(C)};
		final Caps loaded=Caps.load(prefix+".bounds.tsv", "fixture", 3);
		loaded.validateKeys(keys);
		check(loaded.chainCeiling==24 && loaded.chainCeilingRule.equals(NcrnaOverlapWindowBuilder.CEILING_MAX)
			&& loaded.unmeasuredKeys==1 && loaded.leftMargin(keys[2])==-1
			&& loaded.constructionPopulation.equals(NcrnaOverlapWindowBuilder.POPULATION_ORIGINALS),
			"Four original lengths max41 declare ceiling24 while the unseen key stays explicitly unavailable");
		final String scored=root+"/scored_bounds";
		ReportNcrnaSeedPositions.main(new String[]{"mode=bounds", "metadata=t", "population=scored_reference_intervals", "in="+genes,
			"seeds="+seeds, "member=fixture", "expectedgenes=4", "expectedkeys=3", "orientation=sense", "out="+scored, "t=1"});
		check(Caps.load(scored+".bounds.tsv", "fixture", 3).constructionPopulation.equals(NcrnaOverlapWindowBuilder.POPULATION_SCORED),
			"Explicit scored-interval construction retains its exact in-corpus declaration, not an inferred original-gene claim");
		final String regenerated=root+"/regenerated_legacy";
		ReportNcrnaSeedPositions.main(new String[]{"mode=bounds", "in="+genes, "seeds="+seeds,
			"member=fixture", "expectedgenes=4", "expectedkeys=3", "orientation=sense", "out="+regenerated, "t=1"});
		final Caps legacyProduced=Caps.load(regenerated+".bounds.tsv", "fixture", 3);
		CobalaminCallerResources.requireLegacyBounds(legacyProduced);
		check(legacyProduced.leftMargin(encode(A))==24 && legacyProduced.rightMargin(encode(A))==24
			&& legacyProduced.leftMargin(encode(B))==24 && legacyProduced.rightMargin(encode(B))==12,
			"Default p99+10 on four literal genes retains the original24/24and24/12margins after the full-member cap");
		check(legacyProduced.chainCeiling==0 && legacyProduced.constructionPopulation==null,
			"The actual default bounds producer still regenerates a legacy-policy-consumable table");
		try(Input in=new Input(regenerated+".bounds.tsv")){
			final LineParser1 fields=new LineParser1('\t');
			int rows=0;
			for(byte[] line; (line=in.reader.nextLine())!=null; ){
				fields.set(line);
				check(fields.terms()==8 && line[0]!='#', "Legacy regeneration retains the old eight-column header/key rows only");
				rows++;
			}
			check(rows==4, "Legacy output contains one header and every one of the three fixed-axis keys");
		}
		for(final long key:keys){
			check(legacyProduced.leftMargin(key)==loaded.leftMargin(key) && legacyProduced.rightMargin(key)==loaded.rightMargin(key),
				"Opting into self-contained metadata never changes the actual per-key end bounds");
		}
		final String legacy=root+"/legacy_bounds.tsv";
		final String rows="fixture\t"+A+"\t2\t1\t10\t20\t47\t"+NcrnaOverlapWindowBuilder.CAP_RULE+'\n'
			+"fixture\t"+B+"\t0\t0\tNA\tNA\tNA\tNO_OCCURRENCES\n";
		write(legacy, NcrnaOverlapWindowBuilder.CAP_HEADER+'\n'+rows);
		final Caps old=Caps.load(legacy, "fixture", 2);
		check(old.chainCeiling==0 && old.chainCeilingRule==null && old.unmeasuredKeys==1,
			"Legacy cobalamin tables retain absent metadata and their external policy contract");
		CobalaminCallerResources.requireLegacyBounds(old);
		expect(()->CobalaminCallerResources.requireLegacyBounds(loaded));
		final String bad=root+"/duplicate_metadata.tsv";
		final String metadata=NcrnaOverlapWindowBuilder.CEILING_TAG+"\t24\t"+NcrnaOverlapWindowBuilder.CEILING_MAX+'\n';
		write(bad, NcrnaOverlapWindowBuilder.CAP_HEADER+'\n'+metadata+metadata+rows);
		expect(()->Caps.load(bad, "fixture", 2));
		final String badPopulation=root+"/bad_population.tsv";
		write(badPopulation, NcrnaOverlapWindowBuilder.CAP_HEADER+'\n'+metadata
			+NcrnaOverlapWindowBuilder.POPULATION_TAG+"\tunverified guess\n"+rows);
		expect(()->Caps.load(badPopulation, "fixture", 2));
		expect(()->old.validateKeys(new long[]{encode(A), encode(C)}));
	}
	/** Literal tail distances distinguish both requested wider rules and exercise actual generic scheduling. */
	private static void parameterBounds(final String root){
		final int[] distances=new int[1000];
		Arrays.fill(distances, 2);
		Arrays.fill(distances, 990, 999, 30);
		distances[999]=80;
		check(ReportNcrnaSeedPositions.boundedMargin(distances, 1000, 200)==12
			&& ReportNcrnaSeedPositions.boundedMargin(distances, 1000, 200, 99, 20)==22
			&& ReportNcrnaSeedPositions.boundedMargin(distances, 1000, 200, 99.9, 10)==40,
			"With990distances2,nine30andone80, nearest-rankp99+10/p99+20/p99.9+10 must be12/22/40bases");
		check(ReportNcrnaSeedPositions.boundedMargin(new int[]{30}, 1, 50, 99.9, 10)==33,
			"Added padding still caps at the full member maximum50minus17, not a percentile-specific maximum");
		expect(()->NcrnaOverlapWindowBuilder.capRule(99, -1));
		expect(()->NcrnaOverlapWindowBuilder.capRule(Double.NaN, 10));
		expect(()->NcrnaOverlapWindowBuilder.capRule(99.95, 10));
		final String genes=root+"/tail_genes.fa", seeds=root+"/tail_seeds.fa";
		final byte[] sequence=new byte[200];
		final byte[] observed=bytes(A);
		final ByteBuilder buffer=new ByteBuilder();
		try(Output out=new Output(genes)){
			for(int i=0; i<1000; i++){
				Arrays.fill(sequence, (byte)'N');
				System.arraycopy(observed, 0, sequence, distances[i], 17);
				out.writer.print(buffer.clear().append(">tail").append(i).nl().append(sequence).nl());
			}
		}
		write(seeds, ">observed\n"+A+"\n>unseen\n"+C+'\n');
		final String[] percentiles={"99", "99", "99.9"};
		final int[] padding={10, 20, 10}, expected={12, 22, 40};
		final String[] rules={"P99_PLUS10_MEMBER_MAX_MINUS17", "P99_PLUS20_MEMBER_MAX_MINUS17", "P99.9_PLUS10_MEMBER_MAX_MINUS17"};
		final LongHashSet keys=new LongHashSet(4);
		keys.add(encode(A));
		keys.add(encode(C));
		final byte[][] library={bytes(A)};
		final String[] names={"fixture_model"};
		final NcrnaFamily family=new NcrnaFamily("fixture", library, null, names, keys, 17, 40, 83,
			7, 1, false, 0f, 0f, 0f, 1, 0f, 20f, .75f, .65f);
		for(int i=0; i<3; i++){
			final String prefix=root+"/tail_bounds_"+i;
			final ArrayList<String> arguments=new ArrayList<String>(Arrays.asList("mode=bounds", "metadata=t", "in="+genes, "seeds="+seeds,
				"member=fixture", "expectedgenes=1000", "expectedkeys=2", "orientation=sense",
				"out="+prefix, "t=1"));
			if(i>0){
				arguments.add("boundpercentile="+percentiles[i]);
				arguments.add("boundpadding="+padding[i]);
			}
			ReportNcrnaSeedPositions.main(arguments.toArray(new String[arguments.size()]));
			verifyBoundRule(prefix+".bounds.tsv", expected[i], rules[i]);
			final Caps caps=Caps.loadGeneric(prefix+".bounds.tsv", "fixture", 2);
			check(caps.chainCeiling==183 && caps.unmeasuredKeys==1 && caps.leftMargin(encode(A))==expected[i]
				&& caps.rightMargin(encode(A))==183 && caps.leftMargin(encode(C))==-1,
				"Actual producer/consumer preserves the fixed axis, unavailable key, full-member cap and unchanged chain ceiling");
			family.keyBounds=caps;
			final NcrnaScavenger worker=new NcrnaScavenger(library, null, names, keys, 17, 40, 83,
				7, 1, false, 0f, 0f, 0f, 1, 0f, 20f, .75f, .65f);
			GeneCaller.applyNcrnaFamilyControls(worker, family);
			final ArrayList<int[]> windows=worker.prepareCandidateWindows(new int[]{200}, new long[]{encode(A)}, 500);
			check(windows.size()==1 && windows.get(0)[0]==192-expected[i] && windows.get(0)[1]==391,
				"The actual generic worker must consume the selected12/22/40base left margin, preserving its wider original interval");
			if(i>0){expect(()->Caps.load(prefix+".bounds.tsv", "fixture", 2));}
		}
		final String mixed=root+"/mixed_bound_rules.tsv";
		final String header=NcrnaOverlapWindowBuilder.CAP_HEADER+'\n'+NcrnaOverlapWindowBuilder.CEILING_TAG+"\t183\t"
			+NcrnaOverlapWindowBuilder.CEILING_P999+'\n';
		write(mixed, header+"fixture\t"+A+"\t1\t1\t12\t20\t49\tP99_PLUS20_MEMBER_MAX_MINUS17\n"
			+"fixture\t"+B+"\t1\t1\t12\t20\t49\tP99.9_PLUS10_MEMBER_MAX_MINUS17\n");
		expect(()->Caps.loadGeneric(mixed, "fixture", 2));
		final String legacy=root+"/old_diagnostic_bounds.tsv";
		write(legacy, header+"fixture\t"+A+"\t1\t1\t10\t20\t47\tNEAREST_RANK_P999\n");
		expect(()->Caps.loadGeneric(legacy, "fixture", 1));
		expect(()->Caps.loadGeneric(root+"/regenerated_legacy.bounds.tsv", "fixture", 3));
	}
	/** Literal distinct/repeated keys distinguish the actual worker's union and intersection rules. */
	private static void combiningBounds(final String root){
		final String prefix=root+"/intersection_produced";
		ReportNcrnaSeedPositions.main(new String[]{"mode=bounds", "metadata=t", "boundcombine=intersection",
			"in="+root+"/tail_genes.fa", "seeds="+root+"/tail_seeds.fa", "member=fixture", "expectedgenes=1000",
			"expectedkeys=2", "orientation=sense", "out="+prefix, "t=1"});
		final Caps produced=Caps.loadGeneric(prefix+".bounds.tsv", "fixture", 2);
		check(!produced.union && produced.leftMargin(encode(A))==12 && produced.chainCeiling==183,
			"The actual producer must declare intersection without changing default percentile, padding or protection");
		expect(()->Caps.load(prefix+".bounds.tsv", "fixture", 2));
		final LongHashSet keys=new LongHashSet(4);
		keys.add(encode(A));
		keys.add(encode(B));
		final byte[][] library={bytes(A)};
		final String[] names={"fixture_model"};
		final byte[] sequence=new byte[500];
		Arrays.fill(sequence, (byte)'N');
		System.arraycopy(bytes(A), 0, sequence, 92, 17);
		System.arraycopy(bytes(B), 0, sequence, 112, 17);
		for(final String combination:new String[]{"UNION", "INTERSECTION"}){
			final String path=root+"/combination_"+combination+".tsv";
			final String metadata=NcrnaOverlapWindowBuilder.CAP_HEADER+'\n'+NcrnaOverlapWindowBuilder.CEILING_TAG
				+"\t200\t"+NcrnaOverlapWindowBuilder.CEILING_MAX+'\n'
				+NcrnaOverlapWindowBuilder.COMBINATION_TAG+'\t'+combination+'\n';
			write(path, metadata+"fixture\t"+A+"\t1\t1\t40\t80\t137\tP99_PLUS10_MEMBER_MAX_MINUS17\n"
				+"fixture\t"+B+"\t1\t1\t40\t80\t137\tP99_PLUS10_MEMBER_MAX_MINUS17\n");
			final Caps caps=Caps.loadGeneric(path, "fixture", 2);
			final boolean union=combination.equals("UNION");
			check(caps.union==union, "The loaded combining rule must equal its literal metadata value");
			final NcrnaFamily family=new NcrnaFamily("fixture", library, null, names, keys, 17, 40, 83,
				7, 1, false, 0f, 0f, 0f, 1, 0f, 20f, .75f, .65f);
			family.keyBounds=caps;
			final NcrnaScavenger worker=new NcrnaScavenger(library, null, names, keys, 17, 40, 83,
				7, 1, false, 0f, 0f, 0f, 1, 0f, 20f, .75f, .65f);
			GeneCaller.applyNcrnaFamilyControls(worker, family);
			final int start=union ? 52 : 72, stop=union ? 208 : 188;
			final ArrayList<int[]> windows=worker.prepareCandidateWindows(new int[]{100, 120}, new long[]{encode(A), encode(B)}, 500);
			check(windows.size()==1 && windows.get(0)[0]==start && windows.get(0)[1]==stop && windows.get(0)[3]==200,
				"Distinct keys propose52..188and72..208: union52..208versusintersection72..188, with protection200");
			final ArrayList<int[]> repeats=worker.prepareCandidateWindows(new int[]{100, 120}, new long[]{encode(A), encode(A)}, 500);
			check(repeats.size()==1 && repeats.get(0)[0]==52 && repeats.get(0)[1]==208,
				"Repeated occurrences of one key retain their union envelope under either distinct-key combining rule");
			final TrnaCaller trna=new TrnaCaller(new GeneModel(true));
			trna.setKeyBounds(caps);
			final ArrayList<int[]> independent=trna.prepareCandidateWindows(sequence, new int[]{100, 120});
			check(independent.size()==1 && independent.get(0)[0]==start && independent.get(0)[1]==stop,
				"The independent scheduling seam honors the same generic table without changing its verifier");
			final String duplicate=root+"/duplicate_"+combination+".tsv";
			write(duplicate, metadata+NcrnaOverlapWindowBuilder.COMBINATION_TAG+"\tUNION\n");
			expect(()->Caps.loadGeneric(duplicate, "fixture", 2));
			expect(()->Caps.load(path, "fixture", 2));
		}
		final String malformed=root+"/bad_combination.tsv";
		write(malformed, NcrnaOverlapWindowBuilder.CAP_HEADER+'\n'+NcrnaOverlapWindowBuilder.COMBINATION_TAG+"\tintersection\n");
		expect(()->Caps.loadGeneric(malformed, "fixture", 2));
	}
	/** Direct row oracles keep a mutually wrong producer/parser from agreeing on a mislabeled estimator. */
	private static void verifyBoundRule(final String path, final int expectedLeft, final String rule){
		final LineParser1 fields=new LineParser1('\t');
		int rows=0;
		try(Input in=new Input(path)){
			in.reader.nextLine();
			for(byte[] line; (line=in.reader.nextLine())!=null; ){
				if(line[0]=='#'){continue;}
				fields.set(line);
				check(fields.terms()==8, "Parameterized bounds preserve the eight-column serving schema");
				if(fields.termEquals(A, 1)){
					check(fields.parseInt(2)==1000 && fields.parseInt(3)==1000 && fields.parseInt(4)==expectedLeft
						&& fields.parseInt(5)==183 && fields.termEquals(rule, 7), "The declared rank/padding must match literal observed end distances");
				}else{
					check(fields.termEquals(C, 1) && fields.parseInt(2)==0 && fields.termEquals("NO_OCCURRENCES", 7),
						"Parameter changes cannot invent observations for an unobserved forward key");
				}
				rows++;
			}
		}
		check(rows==2, "Both fixed-axis keys survive every parameterized table");
	}

	/** Literal tandem windows exercise the production shared builder and independent tRNA route. */
	private static void windowGeometry(){
		final long[] keys={encode(A), encode(B), encode(C)};
		final int[] centers={100, 130, 160};
		final Caps caps=new Caps(keys, new int[]{10, 10, 10}, new int[]{20, 20, 20}, 100, NcrnaOverlapWindowBuilder.CEILING_MAX);
		final NcrnaOverlapWindowBuilder builder=new NcrnaOverlapWindowBuilder();
		final ArrayList<int[]> windows=builder.build(centers, keys, 400, 83, true, 1, 0, caps, 100, true);
		check(windows.size()==2 && windows.get(0)[0]==82 && windows.get(0)[1]==158
			&& windows.get(1)[0]==142 && windows.get(1)[1]==188,
			"A three-hit overlap chain splits into77and47base windows before swallowing neighboring copies");
		final ArrayList<int[]> intersections=builder.build(centers, keys, 400, 83, false, 1, 0, caps, 100, true);
		check(intersections.size()==2 && intersections.get(0)[0]==112 && intersections.get(0)[1]==128
			&& intersections.get(1)[0]==142 && intersections.get(1)[1]==188
			&& intersections.get(0)[3]==100 && intersections.get(1)[3]==100,
			"Intersection must use the same protected raw-chain split before tightening distinct-key envelopes");
		final byte[] sequence=new byte[400];
		Arrays.fill(sequence, (byte)'N');
		final String[] words={A, B, C};
		for(int i=0; i<words.length; i++){System.arraycopy(bytes(words[i]), 0, sequence, centers[i]-8, 17);}
		check(Arrays.equals(NcrnaOverlapWindowBuilder.hitKeys(sequence, centers, 17), keys),
			"Reconstructed supplied hits retain the exact forward scanner axis without another genome scan");
		final TrnaCaller trna=new TrnaCaller(new GeneModel(true));
		final ArrayList<int[]> legacy=trna.prepareCandidateWindows(sequence, new int[]{100});
		check(legacy.size()==1 && legacy.get(0)[0]==17 && legacy.get(0)[1]==200,
			"No-bounds tRNA retains its old center83padding plus17 inclusive convention");
		trna.setKeyBounds(caps);
		final ArrayList<int[]> typed=trna.prepareCandidateWindows(sequence, centers);
		check(typed.size()==windows.size() && Arrays.equals(typed.get(0), windows.get(0))
			&& Arrays.equals(typed.get(1), windows.get(1)), "Independent tRNA scheduling shares the actual protected geometry");
		final ArrayList<int[]> claims=new ArrayList<int[]>();
		claims.add(new int[]{135, 145});
		final ArrayList<int[]> wide=new ArrayList<int[]>();
		wide.add(new int[]{80, 220, 1, 200});
		final ArrayList<int[]> split=TrnaCaller.subtractClaimed(wide, claims);
		check(split.size()==2 && split.get(0)[0]==80 && split.get(0)[1]==134
			&& split.get(1)[0]==146 && split.get(1)[1]==220
			&& split.get(0).length==4 && split.get(0)[3]==200 && split.get(1)[3]==200,
			"Both claim remainders retain protection, with no verifier or claim-length change");
		final Caps unavailable=new Caps(new long[]{keys[0]}, new int[]{-1}, new int[]{-1}, 60, NcrnaOverlapWindowBuilder.CEILING_MAX);
		final ArrayList<int[]> fallback=builder.build(new int[]{100}, new long[]{keys[0]}, 400, 83, true, 1, 0, unavailable, 60, true);
		check(fallback.size()==1 && fallback.get(0)[0]==40 && fallback.get(0)[1]==177 && builder.unavailableSeeds==1,
			"An explicitly unmeasured generic key keeps old geometry with capped padding");
		check(builder.build(new int[]{100}, new long[]{keys[0]}, 400, 83, true, 1, 0, unavailable, 60).isEmpty(),
			"Existing cobalamin skip-unmeasured behavior stays unchanged");
		final Caps giant=new Caps(new long[]{keys[0]}, new int[]{100}, new int[]{100}, 60, NcrnaOverlapWindowBuilder.CEILING_MAX);
		final ArrayList<int[]> original=builder.build(new int[]{150}, new long[]{keys[0]}, 400, 83, true, 1, 0, giant, 60, true);
		check(original.size()==1 && original.get(0)[0]==42 && original.get(0)[1]==258,
			"An original217base proposal wider than the ceiling is preserved, not clipped to60");
	}
	/** A sibling-only hit cannot open a window; a genuinely unseen key still can. */
	private static void memberAttribution(){
		final long[] keys={encode(A), encode(B), encode(C), encode("GGGGGGGGGGGGGGGGG")};
		final LongHashSet seeds=new LongHashSet(8);
		for(final long key:keys){seeds.add(key);}
		final byte[][] library={bytes(A+"NNN"+B)};
		final String[] names={"fixture_model"};
		final NcrnaFamily small=new NcrnaFamily("srp_small", library, null, names, seeds, 17, 40, 83,
			7, 1, false, 0f, 0f, 0f, 1, 0f, 20f, .75f, .65f);
		final NcrnaFamily large=new NcrnaFamily("srp_large", library, null, names, seeds, 17, 40, 83,
			7, 1, false, 0f, 0f, 0f, 1, 0f, 20f, .75f, .65f);
		final NcrnaFamily unrelated=new NcrnaFamily("rnasep", library, null, names, seeds, 17, 40, 83,
			7, 1, false, 0f, 0f, 0f, 1, 0f, 20f, .75f, .65f);
		small.keyBounds=new Caps(keys, new int[]{10, -1, -1, 15}, new int[]{20, -1, -1, 25}, 60, NcrnaOverlapWindowBuilder.CEILING_MAX);
		large.keyBounds=new Caps(keys, new int[]{-1, 20, -1, 25}, new int[]{-1, 30, -1, 35}, 60, NcrnaOverlapWindowBuilder.CEILING_MAX);
		unrelated.keyBounds=new Caps(keys, new int[]{10, 10, 10, 10}, new int[]{20, 20, 20, 20}, 60, NcrnaOverlapWindowBuilder.CEILING_MAX);
		final Caps originalSmall=small.keyBounds;
		final Caps originalLarge=large.keyBounds;
		final ArrayList<NcrnaFamily> families=new ArrayList<NcrnaFamily>();
		families.add(small);
		families.add(large);
		families.add(unrelated);
		NcrnaFamilyLoader.attributeMemberKeys(families);
		check(small.keyBounds.siblingOnlyKeys==1 && large.keyBounds.siblingOnlyKeys==1
			&& unrelated.keyBounds.siblingOnlyKeys==0,
			"Only observations in the same registered group can exclude an unobserved member key");
		check(originalSmall.siblingOnlyKeys==0 && originalLarge.siblingOnlyKeys==0,
			"Attribution publishes new immutable masks without changing the original loaded caps");
		final NcrnaScavenger worker=new NcrnaScavenger(library, null, names, seeds, 17, 40, 83,
			7, 1, false, 0f, 0f, 0f, 1, 0f, 20f, .75f, .65f);
		GeneCaller.applyNcrnaFamilyControls(worker, small);
		final ArrayList<int[]> windows=worker.prepareCandidateWindows(new int[]{100, 220, 340, 460}, keys, 600);
		check(windows.size()==3 && windows.get(0)[0]==82 && windows.get(0)[1]==128
			&& windows.get(1)[0]==280 && windows.get(1)[1]==417
			&& windows.get(2)[0]==437 && windows.get(2)[1]==493,
			"Own and both-observed keys keep this member's margins, sibling-only opens nothing, and neither-observed keeps capped fallback");
		final NcrnaOverlapWindowBuilder builder=new NcrnaOverlapWindowBuilder();
		check(builder.build(new int[]{220}, new long[]{keys[1]}, 500, 83, true, 1, 0,
			originalSmall, 60, true).size()==1,
			"Standalone unobserved caps retain their previous fallback when no sibling evidence was bound");
		check(builder.build(new int[]{340}, new long[]{keys[2]}, 500, 83, true, 1, 0,
			small.keyBounds, 60, false).isEmpty(),
			"Legacy skip-unobserved behavior is unchanged by the generic attribution mask");
	}
	private static void write(final String path, final String text){
		try(Output out=new Output(path)){out.writer.print(bytes(text));}
	}
	/** Exact per-worker counters cannot lose updates to another worker's calls. */
	private static void parallelTrnaCounts(){
		final byte[] sequence=bytes("AGCTTGACCGTATGGCATTCGAGTACCGATGCATGGCTAGTCGATTCGACCTGATCGTAGCTAGGCTAACGTAGCTTCAGGCTA");
		final byte[][] library={sequence};
		final LongHashSet prior=ProkObject.trnaKmers;
		ProkObject.trnaKmers=null;//Exercise the verifier directly, independent of discovery resources.
		try{
			final TrnaCaller baseline=new TrnaCaller(new GeneModel(true), library);
			baseline.alignWindow("fixture", sequence, 0, 0, sequence.length-1, new ArrayList<Orf>());
			check(baseline.workerAlignmentCount()==1, "One exact model in a short window makes one physical verifier request");
			final TrnaCaller[] callers={new TrnaCaller(new GeneModel(true), library), new TrnaCaller(new GeneModel(true), library)};
			final Throwable[] failures=new Throwable[2];
			final Thread[] workers=new Thread[2];
			final int rounds=8;
			for(int i=0; i<workers.length; i++){
				final int worker=i;
				workers[i]=new Thread(()->{
					try{
						final ArrayList<Orf> output=new ArrayList<Orf>();
						for(int j=0; j<rounds; j++){
							output.clear();
							callers[worker].alignWindow("fixture", sequence, 0, 0, sequence.length-1, output);
						}
					}catch(Throwable failure){failures[worker]=failure;}
				});
				workers[i].start();
			}
			for(final Thread worker:workers){
				try{worker.join();}
				catch(InterruptedException failure){
					Thread.currentThread().interrupt();
					throw new AssertionError("All verifier workers must finish before counting their work", failure);
				}
			}
			for(int i=0; i<workers.length; i++){
				check(failures[i]==null, "Parallel verifier fixture failed: "+failures[i]);
				check(callers[i].workerAlignmentCount()==rounds && callers[i].windowCount()==0,
					"Each worker owns eight verifier calls; direct verification is not window scheduling");
			}
			check(callers[0].workerAlignmentCount()+callers[1].workerAlignmentCount()==2L*rounds,
				"The report merges independent primitive totals after joining, without a shared atomic");
		}finally{ProkObject.trnaKmers=prior;}
	}

	/** Actual ordinary verification consumes family-bound resources after subtracting a prior claim. */
	private static void ordinaryWorker(){
		final byte[] bases=new byte[400];
		Arrays.fill(bases, (byte)'G');
		final int[] centers={100, 130, 160};
		final String[] words={A, B, C};
		final long[] keys={encode(A), encode(B), encode(C)};
		final LongHashSet seeds=new LongHashSet(8);
		for(int i=0; i<words.length; i++){
			System.arraycopy(bytes(words[i]), 0, bases, centers[i]-8, 17);
			seeds.add(keys[i]);
		}
		// The surviving first query is exactly this58base model; the second need not pass.
		final byte[][] library={Arrays.copyOfRange(bases, 101, 159)};
		final String[] names={"fixture_model"};
		final NcrnaFamily family=new NcrnaFamily("fixture", library, null, names, seeds, 17, 40, 83,
			7, 1, false, 0f, 0f, 0f, 1, 0f, 20f, .75f, .65f);
		family.keyBounds=new Caps(keys, new int[]{10, 10, 10}, new int[]{20, 20, 20}, 100, NcrnaOverlapWindowBuilder.CEILING_MAX);
		family.scavengePass2=false;
		family.trimAlignmentExtent=false;
		family.modelPatience=-1;
		final ArrayList<int[]> scheduled=new ArrayList<int[]>();
		final NcrnaScavenger worker=new NcrnaScavenger(library, null, names, seeds, 17, 40, 83,
			7, 1, false, 0f, 0f, 0f, 1, 0f, 20f, .75f, .65f){
			@Override
			protected void scheduledWindow(final int[] window){scheduled.add(window.clone());}
		};
		final boolean refresh=CallGenes.NCRNA_REFRESH_CLAIMED_WINDOWS;
		CallGenes.NCRNA_REFRESH_CLAIMED_WINDOWS=false;
		try{GeneCaller.applyNcrnaFamilyControls(worker, family);}
		finally{CallGenes.NCRNA_REFRESH_CLAIMED_WINDOWS=refresh;}
		final ArrayList<int[]> claims=new ArrayList<int[]>();
		claims.add(new int[]{0, 100});
		final ArrayList<Orf> calls=worker.scavenge("fixture", bases, 0, claims, centers);
		check(!calls.isEmpty() && worker.alignmentCount()>0, "The ordinary family-controlled worker must perform positive real verification");
		check(scheduled.size()==2 && scheduled.get(0)[0]==101 && scheduled.get(0)[1]==158
			&& scheduled.get(1)[0]==142 && scheduled.get(1)[1]==188,
			"The scheduled queries are58and47bases after the initial claim, not their raw77and47base proposals");
		check(scheduled.get(0).length==4 && scheduled.get(0)[3]==100 && scheduled.get(1)[3]==100,
			"The family hook and ordinary claim subtraction retain the declared neighbor protection");
		check(worker.windowCount()==2 && worker.windowBases()==105,
			"Physical work is two post-claim windows totalling105bases, independent of how many models align or accept");
	}
	/** Resolver-only fixtures need no biological payload or extra program copy. */
	private static void resourceBinding(final String root){
		final String directory=root+"/standalone";
		check(new File(directory).mkdir(), "Resource fixture directory must be fresh");
		final String header="#ncrna_family\t1\tsam1\nrole\tmodel\tend\tpath\tformat\tbytes\tsha80\n";
		write(directory+"/manifest.tsv", header);
		final Overrides options=new Overrides();
		options.root=directory;
		check(new FamilyResourceResolver().optionalKeyBounds("sam1", options)==null,
			"An absent optional table does not invent a serving policy");
		write(directory+"/key_bounds.tsv", "fixture payload\n");
		check(new FamilyResourceResolver().optionalKeyBounds("sam1", options).path.equals(directory+"/key_bounds.tsv"),
			"The canonical standalone filename enables the resolver binding");
		write(directory+"/key_bounds.tsv.gz", "fixture payload\n");
		expect(()->new FamilyResourceResolver().optionalKeyBounds("sam1", options));
		final Overrides declared=new Overrides();
		declared.root=root+"/declared";
		check(new File(declared.root).mkdir(), "Explicit resource fixture must be fresh");
		write(declared.root+"/declared.tsv", "fixture payload\n");
		write(declared.root+"/key_bounds.tsv", "fixture payload\n");
		write(declared.root+"/key_bounds.tsv.gz", "fixture payload\n");
		write(declared.root+"/manifest.tsv", header+"key_bounds\t.\t.\tdeclared.tsv\ttsv\t-\t-\n");
		check(new FamilyResourceResolver().optionalKeyBounds("sam1", declared).path.equals(declared.root+"/declared.tsv"),
			"An explicit manifest selects one table even when both canonical spellings exist");
		final Overrides missing=new Overrides();
		missing.root=root+"/missing";
		check(new File(missing.root).mkdir(), "Missing-payload fixture must be fresh");
		write(missing.root+"/manifest.tsv", header+"key_bounds\t.\t.\tabsent.tsv\ttsv\t-\t-\n");
		expect(()->new FamilyResourceResolver().optionalKeyBounds("sam1", missing));
		final NcrnaFamilyLoader.Options disabled=new NcrnaFamilyLoader.Options();
		disabled.enabled=false;
		disabled.resources.root=root+"/must_not_be_opened";
		check(NcrnaFamilyLoader.load("sam1", disabled, new FamilyResourceResolver())==null,
			"Disabled families open neither optional tables nor a missing manifest");
		final String group=root+"/group";
		check(new File(group+"/small").mkdirs() && new File(group+"/large").mkdir(), "Grouped member fixtures are separate data directories");
		write(group+"/manifest.tsv", "#ncrna_family\t2\tSRP\nmember\trole\tmodel\tend\tpath\tformat\tbytes\tsha80\n"
			+".\tseeds\t.\t.\tseeds.fa\tfasta\t-\t-\n"
			+"small\tconsensus\t.\t.\tconsensus.fa\tfasta\t-\t-\nsmall\tmodels\t.\t.\tmodels.hbm\thbm\t-\t-\n"
			+"large\tconsensus\t.\t.\tconsensus.fa\tfasta\t-\t-\nlarge\tmodels\t.\t.\tmodels.hbm\thbm\t-\t-\n");
		write(group+"/small/key_bounds.tsv", "fixture payload\n");
		final Overrides grouped=new Overrides();
		grouped.root=group;
		final FamilyResourceResolver resolver=new FamilyResourceResolver();
		check(resolver.optionalKeyBounds("srp_small", grouped).path.equals(group+"/small/key_bounds.tsv")
			&& resolver.optionalKeyBounds("srp_large", grouped)==null,
			"Shared discovery seeds do not silently share one member's bounds with another");
	}
	private static void expect(final Runnable action){
		try{action.run();}
		catch(IllegalArgumentException expected){return;}
		throw new AssertionError("Malformed metadata or a changed seed axis must fail before workers start");
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
