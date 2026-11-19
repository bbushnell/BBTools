package test.prok;

import java.io.File;
import java.util.Arrays;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaLearningContext;
import prok.NcrnaLearningPopulation;
import prok.NcrnaPositionalTables;
import prok.NcrnaSupportIO;
import prok.NcrnaTextIO.Input;
import prok.RrnaPositionalKmerTable;
import shared.Shared;
import stream.Read;
import structures.ByteBuilder;

/** Literal all-radius count oracle on the shared21-original population fixture.
 * Source exclusions, missing geometry, missing GC and ambiguity have distinct
 * known dispositions; an unavailable feature never removes its original ID.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaPositionalTablesTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaPositionalTablesTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown positional-table fixture option: "+arg);
			}
			check(parser.out1!=null && parser.in1==null && parser.in2==null && parser.maxReads<0 && parser.samplerate==1
				&& !parser.overwrite && !parser.append && new File(parser.out1).mkdir(), "Supply out=fresh_literal_fixture_directory");
			run(new File(parser.out1));
			pp.outstream.println("NCRNA_POSITIONAL_TABLES_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void run(final File root){
		final String roster=new File(root, "roster.tsv").getPath(), loci=new File(root, "loci.tsv").getPath();
		final String graph=new File(root, "graph.tsv").getPath(), exclude=new File(root, "exclude.txt").getPath();
		final String fasta=new File(root, "contexts.fa").getPath();
		NcrnaLearningPopulationTest.write(roster, NcrnaLearningPopulationTest.roster());
		NcrnaLearningPopulationTest.write(loci, NcrnaLearningPopulationTest.loci());
		NcrnaLearningPopulationTest.write(graph, NcrnaLearningPopulationTest.graph());
		NcrnaLearningPopulationTest.write(exclude, "gene0\ngene2\n");
		NcrnaLearningPopulationTest.write(fasta, corpus());
		final NcrnaLearningPopulation population=NcrnaLearningPopulation.loadForRefit(NcrnaSupportIO.read(graph),
			new String[]{"common", "rare"}, roster, loci, exclude);
		final boolean trim=Shared.TRIM_READ_DESCRIPTION, upper=Read.TO_UPPER_CASE;
		final NcrnaPositionalTables tables=NcrnaPositionalTables.build(population, fasta, "RF00174", new int[]{6, 10, 16}, 1.0/262144);
		final NcrnaPositionalTables defaults=NcrnaPositionalTables.build(population, fasta, "RF00174", new int[]{6, 10, 16},
			9, 0, 9, 0, 1);
		check(Shared.TRIM_READ_DESCRIPTION==trim && Read.TO_UPPER_CASE==upper, "Table counting must restore its borrowed process-wide reader settings");
		check(tables.records==21 && tables.filtered==3 && tables.learning==18 && tables.counted==16 && tables.noGeometry==2,
			"All21originals remain: two excluded plus one partial;18eligible include two without geometry and16 counted loci");
		check(Arrays.equals(tables.availability, new long[]{18, 1, 1, 1}), "All original contexts retain their four literal availability dispositions");
		check(tables.statusByGene[1]==NcrnaLearningContext.NO_GC && tables.statusByGene[5]==NcrnaLearningContext.NO_FLANK
			&& tables.statusByGene[6]==NcrnaLearningContext.MALFORMED_SPAN && tables.statusByGene[20]==NcrnaLearningContext.AVAILABLE,
			"GC-unavailable and filtered/partial loci must keep their original-axis metadata");
		for(int r=0; r<3; r++){
			final int center=tables.radii[r];
			check(tables.five[r].records==16 && tables.three[r].records==16
				&& tables.five[r].counts[center][0]==15 && tables.three[r].counts[center][0]==16,
				"A9 occurs at15true starts and16true stops at every radius; one start is ambiguous and missing GC still teaches tables");
			check(tables.five[r].valid[center]==15 && tables.five[r].ambiguous[center]==1
				&& tables.three[r].valid[center]==16 && tables.three[r].ambiguous[center]==0,
				"A missing word is explicit ambiguity, not a fabricated neutral anchor observation");
			check(tables.five[r].offset==0 && tables.three[r].offset==-8 && tables.five[r].k==9
				&& tables.three[r].k==9 && tables.five[r].alpha==1.0/262144,
				"All radii use fully-inside K9 anchors and the explicit total-prior-one recipe");
			parity(root, "r"+center+".5prime", tables.five[r], defaults.five[r]);
			parity(root, "r"+center+".3prime", tables.three[r], defaults.three[r]);
		}
		final NcrnaPositionalTables placed=NcrnaPositionalTables.build(population, fasta, "RF00174", new int[]{1},
			3, 2, 5, 3, 1);
		check(placed.five[0].k==5 && placed.five[0].offset==-2 && placed.five[0].alpha==1.0/1024
			&& placed.three[0].k==8 && placed.three[0].offset==-4 && placed.three[0].alpha==1.0/65536,
			"Each end binds its own semantic word and preserves total prior one through division by4^K");
		check(placed.five[0].counts[1][640]==15 && placed.three[0].counts[1][63]==16,
			"Literal true START wordGGAAA(key640) occurs15times; STOP wordAAAAATTT(key63) occurs16times");
		final NcrnaPositionalTables triple=NcrnaPositionalTables.buildTriple(population, fasta, "RF00174", new int[]{10},
			7, 2, 8, 1, 1);
		check(triple.triple() && triple.inputs(0)==69 && triple.fiveProbes[0][0]==triple.five[0]
			&& triple.threeProbes[0][0]==triple.three[0] && triple.counted==16,
			"The optional layout aliases its primary tables and counts each eligible locus once per probe");
		final int[] keys={163840, 174762, 0, 3, 0, 262143};
		final long[] counts={15, 16, 16, 16, 16, 16};
		for(int end=0; end<2; end++){
			for(int probe=0; probe<3; probe++){
				final RrnaPositionalKmerTable.Table table=end==0 ? triple.fiveProbes[0][probe] : triple.threeProbes[0][probe];
				check(table.records==16 && table.counts[10][keys[3*end+probe]]==counts[3*end+probe],
					"Literal primary/left/right START words GGAAAAAAA/GGGGGGGGG/AAAAAAAAA and STOP AAAAAAAAT/AAAAAAAAA/TTTTTTTTT must retain their exact original counts");
				check(table.offset==(probe==0 ? (end==0 ? -2 : -7) : probe==1 ? -9 : 1),
					"Exclusive offsets remain signed sense coordinates at both ends");
			}
		}
		fails(()->NcrnaPositionalTables.wordLength(10, 0));
		fails(()->NcrnaPositionalTables.wordLength(6, -1));
		fails(()->new RrnaPositionalKmerTable.Table("RF00174", "5prime", 10, 0, 1.0/1048576, 1));
	}

	/** The retained per-key-alpha API is the exact original K9 construction route. */
	private static void parity(final File root, final String name, final RrnaPositionalKmerTable.Table before,
			final RrnaPositionalKmerTable.Table after){
		final String first=new File(root, name+".legacy.tsv").getPath(), second=new File(root, name+".default.tsv").getPath();
		before.write(first);
		after.write(second);
		try(Input a=new Input(first); Input b=new Input(second)){
			for(byte[] line; (line=a.reader.nextLine())!=null;){
				check(Arrays.equals(line, b.reader.nextLine()), "Default geometry must preserve every serialized K9 table byte");
			}
			check(b.reader.nextLine()==null, "Default geometry cannot add trailing table rows");
		}
	}
	private static void fails(final Runnable action){
		try{action.run();}catch(IllegalArgumentException expected){return;}
		throw new AssertionError("Invalid single-table geometry must fail before allocating dense counts");
	}

	/** Each admitted core has100bases between20/30-base flanks. */
	private static String corpus(){
		final ByteBuilder result=new ByteBuilder();
		for(int gene=0; gene<21; gene++){
			result.append(">gene").append(gene).append(" source=genome").append(gene==20 ? 0 : gene).append(" lflank=20");
			if(gene!=5){result.append(" rflank=30");}
			if(gene!=1){result.append(" contig_gc=0.5");}
			result.nl();
			final byte[] bases=new byte[gene==6 ? 10 : 150];
			if(gene==6){Arrays.fill(bases, (byte)'A');}
			else{
				Arrays.fill(bases, 0, 20, (byte)'G');
				Arrays.fill(bases, 20, 120, (byte)'A');
				Arrays.fill(bases, 120, 150, (byte)'T');
				if(gene==3){bases[20]='N';}
			}
			result.append(bases).nl();
		}
		return result.toString();
	}
	private static void check(final boolean valid, final String reason){if(!valid){throw new AssertionError(reason);}}
}
