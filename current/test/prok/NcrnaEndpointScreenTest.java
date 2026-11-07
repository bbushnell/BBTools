package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaEndpointScreen;
import prok.NcrnaEndpointScreen.Counts;
import prok.NcrnaEndpointTable;
import prok.NcrnaTextIO.Input;
import prok.TrnaBoundaryFeatures.BoundaryType;
import shared.Shared;
import structures.ByteBuilder;

/** Literal count-ratio, repeated-key LOO, semantic geometry and complete-split controls.
 * The CLI fixture has twenty complete homopolymer loci plus one partial copy;
 * expected eighteen training/two validation loci do not depend on shuffled IDs.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaEndpointScreenTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaEndpointScreenTest.class, false);
		try{
			final Parser parser=new Parser();
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(java.util.Locale.ROOT), b=split.length>1 ? split[1] : null;
				if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown fixture argument: "+arg);}
			}
			check(parser.out1!=null, "Specify out=fresh_directory for the fixture");
			final File root=new File(parser.out1);
			check(root.mkdir(), "Fixture output must be fresh");
			counts();
			population(root);
			System.out.println("NCRNA_ENDPOINT_SCREEN_TEST_PASS native_LOO ratio_smoothing ties semantic_geometry missing_ambiguous split_conservation");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Hand-built words make every raw count and ratio independently calculable. */
	private static void counts(){
		final double[] ratios={1.5, 2, 3, 5, 10};
		for(final BoundaryType type:new BoundaryType[]{BoundaryType.START, BoundaryType.STOP}){
			final boolean stop=type==BoundaryType.STOP;
			final byte[] bases=new byte[36];
			Arrays.fill(bases, (byte)'A');
			// START(inside1,outside2) at8 uses6..8; STOP(inside2,outside1) at27 uses26..28.
			final int word=stop ? 26 : 6, boundary=stop ? 27 : 8;
			System.arraycopy("CGT".getBytes(StandardCharsets.US_ASCII), 0, bases, word, 3);
			final NcrnaEndpointTable table=new NcrnaEndpointTable(type, stop ? 2 : 1, stop ? 1 : 2);
			for(int i=0; i<4; i++){table.add(bases, boundary);}
			final Counts train=new Counts(table, false, 2, ratios), validation=new Counts(table, true, 2, ratios);
			train.observe(bases, 8, 27);
			validation.observe(bases, 8, 27);
			train.validate();
			validation.validate();
			check(train.trueSum==3 && validation.trueSum==4 && train.otherSum==0 && validation.otherSum==0,
				"Only training removes its own unique CGT occurrence; semantic end geometry must hit the hand-placed word");
			check(Arrays.equals(train.correct, new long[]{1, 1, 1, 0, 0})
				&& Arrays.equals(validation.correct, new long[]{1, 1, 1, 1, 0})
				&& Arrays.equals(train.wrong, new long[5]) && train.strictBest==1,
				"Smoothed ratios are exactly4 and5 against a zero-count neighbor; the threshold comparison is inclusive");
			final byte[] shifted=new byte[36];
			Arrays.fill(shifted, (byte)'A');
			System.arraycopy("CGT".getBytes(StandardCharsets.US_ASCII), 0, shifted, word+1, 3);
			final Counts wrong=new Counts(table, true, 2, ratios);
			wrong.observe(shifted, 8, 27);
			wrong.validate();
			check(Arrays.equals(wrong.wrong, new long[]{1, 1, 1, 1, 0}) && wrong.trueSum==0 && wrong.otherSum==4,
				"A neighboring CGT with raw count4 decisively beats the unseen true word by smoothed ratio5");
		}

		final byte[] repeated=new byte[30];
		Arrays.fill(repeated, (byte)'A');
		final NcrnaEndpointTable table=new NcrnaEndpointTable(BoundaryType.START, 3, 0);
		table.add(repeated, 5);
		table.add(repeated, 5);
		final Counts repeatedKey=new Counts(table, false, 2, ratios);
		repeatedKey.observe(repeated, 5, 24);
		repeatedKey.validate();
		check(repeatedKey.trueSum==1 && repeatedKey.otherSum==1 && repeatedKey.tiedBest==1
			&& Arrays.equals(repeatedKey.correct, new long[5]) && Arrays.equals(repeatedKey.wrong, new long[5]),
			"Self subtraction must follow AAA at every matching occurrence, so a repeated-key tie abstains");
		final NcrnaEndpointTable outside=new NcrnaEndpointTable(BoundaryType.START, 1, 2);
		final Counts edge=new Counts(outside, true, 2, ratios);
		edge.observe(repeated, 0, 19);
		edge.validate();
		check(edge.trueMissing==1 && edge.invalidSpan==2 && edge.missingWords==2 && edge.validWords==1 && edge.eligible==0,
			"At start0, two negative offsets leave the sequence and windows at offsets0/1 lack outside bases");
		final byte[] ambiguous=repeated.clone();
		ambiguous[7]='N';
		final Counts ns=new Counts(outside, true, 2, ratios);
		ns.observe(ambiguous, 8, 27);
		ns.validate();
		check(ns.trueAmbiguous==1 && ns.ambiguousWords==3 && ns.validWords==2 && ns.eligible==0,
			"An N in the true word invalidates three adjacent words, without pretending they were valid zero-count competitors");
		final byte[] exactSpan=new byte[16];
		Arrays.fill(exactSpan, (byte)'A');
		final Counts alone=new Counts(table, true, 1, ratios);
		alone.observe(exactSpan, 0, 15);
		alone.validate();
		check(alone.noNeighbors==1 && alone.invalidSpan==2 && alone.validWords==1 && alone.eligible==0,
			"A sixteen-base span without flanks has a true word but no valid displaced gene span");
	}

	/** Exercises real original-population splitting, table construction and native output. */
	private static void population(final File root){
		final String roster=new File(root, "roster.tsv").getPath(), loci=new File(root, "loci.tsv").getPath();
		final String graph=new File(root, "graph.tsv").getPath(), original=new File(root, "original.fa").getPath();
		final String fasta=new File(root, "contexts.fa").getPath();
		NcrnaLearningPopulationTest.write(roster, NcrnaLearningPopulationTest.roster());
		NcrnaLearningPopulationTest.write(loci, NcrnaLearningPopulationTest.loci());
		NcrnaLearningPopulationTest.write(graph, NcrnaLearningPopulationTest.graph());
		NcrnaLearningPopulationTest.write(original, ">common\nAAAAAAAAAAAAAAAAA\n>rare\nCCCCCCCCCCCCCCCCC\n");
		final ByteBuilder b=new ByteBuilder();
		for(int gene=0; gene<21; gene++){
			b.append(">gene").append(gene).append(" source=genome").append(gene==20 ? 0 : gene)
				.append(" lflank=20 rflank=30 contig_gc=0.5\n");
			for(int base=0; base<150; base++){b.append('A');}
			b.nl();
		}
		NcrnaLearningPopulationTest.write(fasta, b.toString());
		final File output=new File(root, "screen");
		NcrnaEndpointScreen.main(new String[]{"in="+fasta, "out="+output, "graph="+graph, "original="+original,
			"roster="+roster, "loci="+loci, "splitseed=42", "ks=3", "radius=1", "ratios=1.0000001,2", "t=2"});
		final LineParser1 row=new LineParser1('\t');
		final boolean[] seen=new boolean[16];
		try(Input input=new Input(new File(output, "summary.tsv").getPath())){
			check(input.reader.nextLine()!=null, "Summary header must exist");
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				row.set(line);
				final boolean validation=row.termEquals("VALIDATION", 4), stop=row.termEquals("STOP", 3);
				final int n=validation ? 2 : 18, outside=row.parseInt(2), key=outside*4+(validation ? 2 : 0)+(stop ? 1 : 0);
				check(row.terms()==30 && row.parseInt(0)==3 && row.parseInt(1)+outside==3 && !seen[key], "Each geometry/end/partition occurs once");
				seen[key]=true;
				check(row.parseInt(6)==(validation ? 2 : 19) && row.parseInt(7)==(validation ? 0 : 1)
					&& row.parseInt(11)==n && row.parseInt(12)==18 && row.parseInt(13)==18,
					"All twenty-one original loci stay visible; the partial is excluded and only eighteen training loci teach");
				check(row.parseInt(16)==n && row.parseInt(20)==0 && row.parseInt(21)==n && row.parseInt(22)==0
					&& row.parseLong(23)==(validation ? 36 : 306) && row.parseLong(24)==(validation ? 36 : 306)
					&& row.parseLong(25)==3*n && row.parseLong(26)==3*n,
					"Training AAA counts are17 everywhere after LOO, validation counts18; all three candidate words tie");
			}
		}
		for(final boolean found:seen){check(found, "All four compositions, both ends and both partitions must be reported");}
		int curves=0;
		try(Input input=new Input(new File(output, "curves.tsv").getPath())){
			input.reader.nextLine();
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				row.set(line);
				check(row.parseDouble(6)==1.0000001 || row.parseDouble(6)==2,
					"A user ratio close to one must retain its value in the report, not round down to one");
				check(row.terms()==13 && row.parseLong(8)==0 && row.parseLong(9)==0 && row.parseLong(10)==row.parseLong(7),
					"Equal counts abstain even at ratio1.0000001; true-position-only self subtraction would create false decisions");
				curves++;
			}
		}
		check(curves==32, "Four geometries times two ends times two partitions times two ratios");
	}

	private static void check(final boolean condition, final String reason){if(!condition){throw new AssertionError(reason);}}
}
