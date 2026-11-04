package test.prok;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import fileIO.ByteStreamWriter;
import prok.NcrnaLearningPopulation;
import prok.NcrnaNumbers;
import prok.NcrnaSupport.Graph;
import prok.NcrnaSupportIO;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;

/** Finite complete-population joins; no biological input or alignment required.
 * Expected group/eligibility counts are independent of the shuffled split.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaLearningPopulationTest{

	public static void main(String[] args) throws Exception{
		assert(args.length==0) : "Population fixture creates all of its own inputs";
		decimalCounts();
		final Path root=Files.createTempDirectory("learning-population-");
		try(AutoCloseable cleanup=()->{
			for(String name:new String[]{"roster.tsv", "loci.tsv", "graph.tsv", "exclude.txt"}){
				Files.deleteIfExists(root.resolve(name));
			}
			Files.delete(root);
		}){
			final String roster=roster(), loci=loci();
			final String rosterFile=root.resolve("roster.tsv").toString();
			final String lociFile=root.resolve("loci.tsv").toString();
			final String graphFile=root.resolve("graph.tsv").toString();
			final String excludedFile=root.resolve("exclude.txt").toString();
			write(rosterFile, roster);
			write(lociFile, loci);
			write(graphFile, graph());
			write(excludedFile, "gene0\ngene2\n");
			final Graph graph=NcrnaSupportIO.read(graphFile);
			final String[] models={"common", "rare"};
			final NcrnaLearningPopulation all=NcrnaLearningPopulation.load(graph, models, rosterFile, lociFile, null, 42);
			final NcrnaLearningPopulation filtered=NcrnaLearningPopulation.load(graph, models, rosterFile, lociFile, excludedFile, 42);
			assert(all.ids.length==21 && all.split().taxids.length==20 && all.split().selected==2)
				: "Twenty positive TaxIDs, one extra gene copy and one empty source yield twenty groups and two holdouts";
			assert(Arrays.equals(all.split().validationByGroup, filtered.split().validationByGroup))
				: "Learning exclusions must not change the whole-population split";
			assert(!all.split().validationByGene[20] && !all.split().validationByGene[0])
				: "The below-cutoff partial gene still protects its rare original-model bin and entire TaxID";
			assert(all.originalModel[20]==1 && all.partial[20] && !all.eligible(20))
				: "Partial loci remain represented with their original winner, but never teach";
			int eligible=0;
			for(int gene=0; gene<filtered.ids.length; gene++){
				if(filtered.eligible(gene)){eligible++;}
				assert(filtered.index("gene"+gene)==gene && filtered.length[gene]==100
					&& filtered.leftFlank[gene]==20 && filtered.rightFlank[gene]==30)
					: "Joined metadata must remain on the graph query axis with sense-oriented flank geometry";
			}
			assert(eligible==18 && filtered.sourceName(20).equals("genome0"))
				: "One partial and two explicit exclusions leave eighteen teaching loci without losing either copy of TaxID1";

			// Reordered ledger still joins by identity, not row position.
			write(lociFile, reversedLoci());
			final NcrnaLearningPopulation reordered=NcrnaLearningPopulation.load(graph, models, rosterFile, lociFile, null, 42);
			assert(Arrays.equals(all.source, reordered.source) && Arrays.equals(all.partial, reordered.partial))
				: "Locus ordering must not change joined source or eligibility columns";
			write(lociFile, loci);
			fails(()->NcrnaLearningPopulation.load(graph, new String[]{"rare", "common"}, rosterFile, lociFile, null, 42));
			fails(()->NcrnaLearningPopulation.load(graph.subset(new boolean[]{true, true}), models, rosterFile, lociFile, null, 42));
			write(lociFile, loci.replace(locus(19), ""));
			fails(()->NcrnaLearningPopulation.load(graph, models, rosterFile, lociFile, null, 42));
			write(lociFile, loci+locus(0));
			fails(()->NcrnaLearningPopulation.load(graph, models, rosterFile, lociFile, null, 42));
			write(lociFile, loci.replace("gene19\t", "unknown\t"));
			fails(()->NcrnaLearningPopulation.load(graph, models, rosterFile, lociFile, null, 42));
			write(lociFile, loci.replace("gene19\tgenome19", "gene19\tgenome0"));
			fails(()->NcrnaLearningPopulation.load(graph, models, rosterFile, lociFile, null, 42));
			write(lociFile, loci.replace("true\tfalse\tSOURCE_PARTIAL", "false\tfalse\tSOURCE_PARTIAL"));
			fails(()->NcrnaLearningPopulation.load(graph, models, rosterFile, lociFile, null, 42));
			write(lociFile, loci.replace("true\tfalse\tSOURCE_PARTIAL", "true\ttrue\tSOURCE_PARTIAL"));
			fails(()->NcrnaLearningPopulation.load(graph, models, rosterFile, lociFile, null, 42));
			write(lociFile, loci.replace("\t20\t30\t1\n", "\t2147483648\t30\t1\n"));
			fails(()->NcrnaLearningPopulation.load(graph, models, rosterFile, lociFile, null, 42));
			write(lociFile, loci);
			write(rosterFile, roster.replace("\t2\t1000\t1\n", "\t1\t1000\t1\n"));
			fails(()->NcrnaLearningPopulation.load(graph, models, rosterFile, lociFile, null, 42));
			write(rosterFile, roster);
			write(excludedFile, "gene2\ngene2\n");
			fails(()->NcrnaLearningPopulation.load(graph, models, rosterFile, lociFile, excludedFile, 42));
			write(excludedFile, "unknown\n");
			fails(()->NcrnaLearningPopulation.load(graph, models, rosterFile, lociFile, excludedFile, 42));
			System.out.println("Learning population passed: complete joins, below-cutoff partial-bin protection, exclusions, reordered loci and twelve invalid inputs");
		}
	}

	/** Literal boundary cases for both text representations, independent of the parser. */
	private static void decimalCounts(){
		final String[] text={"0", "000007", "2147483647"};
		final int[] expected={0, 7, Integer.MAX_VALUE};
		for(int i=0; i<text.length; i++){
			final String framed="x"+text[i]+"y";
			final byte[] bytes=framed.getBytes(StandardCharsets.US_ASCII);
			assert(NcrnaNumbers.unsignedInt(framed, 1, framed.length()-1)==expected[i]
				&& NcrnaNumbers.unsignedInt(bytes, 1, bytes.length-1)==expected[i])
				: "A borrowed decimal range must retain zero, leading zeros and the maximum int: "+text[i];
		}
		for(String invalid:new String[]{"", "-1", "+1", " 1", "1 ", "1.0", "1e2", "2147483648", "18446744073709551616"}){
			final byte[] bytes=invalid.getBytes(StandardCharsets.US_ASCII);
			fails(()->NcrnaNumbers.unsignedInt(invalid, 0, invalid.length()));
			fails(()->NcrnaNumbers.unsignedInt(bytes, 0, bytes.length));
		}
		System.out.println("Strict decimal counts passed: bounded byte/character ranges, zero, leading zeros, maximum int and malformed/overflow rejection");
	}

	static String roster(){
		final ByteBuilder b=new ByteBuilder().append(NcrnaLearningPopulation.ROSTER_HEADER).nl();
		for(int i=0; i<=20; i++){
			b.append(i).tab().append("genome").append(i).tab().append("assembly").append(i)
				.tab().append(i+1).append("\tbacteria\tP\tgff\tfasta\t")
				.append(i==0 ? 2 : i==20 ? 0 : 1).append("\t1000\t1\n");
		}
		return b.toString();
	}

	static String loci(){
		final ByteBuilder b=new ByteBuilder().append(NcrnaLearningPopulation.LOCUS_HEADER).nl();
		for(int i=0; i<=20; i++){b.append(locus(i));}
		return b.toString();
	}

	private static String reversedLoci(){
		final ByteBuilder b=new ByteBuilder().append(NcrnaLearningPopulation.LOCUS_HEADER).nl();
		for(int i=20; i>=0; i--){b.append(locus(i));}
		return b.toString();
	}

	private static String locus(final int gene){
		assert(gene>=0 && gene<=20) : "Fixture gene IDs are zero through twenty";
		final int source=gene==20 ? 0 : gene;
		return new ByteBuilder().append("gene").append(gene).tab().append("genome").append(source)
			.tab().append("assembly").append(source).tab().append(source+1)
			.append("\tbacteria\tP\tcontig\t+\t101-200\t100\t1\tid\t")
			.append(gene==20 ? "true\tfalse\tSOURCE_PARTIAL" : "false\tfalse\tANNOTATION_FULL_COPY")
			.append("\t1000\tfalse\t0\t20\t30\t1\n").toString();
	}

	static String graph(){
		final ByteBuilder b=new ByteBuilder().append("#ncrna_support_graph\t2\n#scope\tORIGINAL\n#aligner\tquantum\n#cutoff_a48\t")
			.appendFloatA48(0.7f).append("\n#models\t2\n#model\t0\tcommon\n#model\t1\trare\n#queries\t21\n")
			.append("query_id\tbest_model_index\tbest_identity_a48\tsupporters\n");
		for(int i=0; i<20; i++){
			b.append("gene").append(i).append("\t0\t").appendFloatA48(0.8f)
				.append("\t0:").appendFloatA48(0.8f).nl();
		}
		b.append("gene20\t1\t").appendFloatA48(0.4f).append("\t.\n");
		return b.toString();
	}

	static void write(final String path, final String text){
		try(Output output=new Output(new ByteStreamWriter(path, true, false, false))){
			output.writer.print(text);
		}
	}

	private static void fails(final Runnable action){
		try{
			action.run();
		}catch(IllegalArgumentException expected){return;}
		throw new AssertionError("Invalid population must fail rather than dropping or relabeling original loci");
	}
}
