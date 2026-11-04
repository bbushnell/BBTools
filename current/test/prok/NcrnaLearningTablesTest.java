package test.prok;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;

import prok.NcrnaLearningContext;
import prok.NcrnaLearningPopulation;
import prok.NcrnaLearningTables;
import prok.NcrnaSupportIO;
import shared.Shared;
import stream.Read;
import structures.ByteBuilder;

/** Complete streamed context join, teaching scope and explicit missing-data cases.
 * Reuses the original-population fixture data, not its expected split algorithm.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaLearningTablesTest{

	public static void main(String[] args) throws Exception{
		assert(args.length==0) : "Context integration fixture has no external corpus";
		final Path root=Files.createTempDirectory("learning-tables-");
		final int buffer=Shared.bufferLen();
		Shared.setBufferLen(2);
		try(AutoCloseable cleanup=()->{
			for(String name:new String[]{"roster.tsv", "loci.tsv", "graph.tsv", "exclude.txt", "contexts.fa"}){
				Files.deleteIfExists(root.resolve(name));
			}
			Files.delete(root);
		}){
			final String roster=root.resolve("roster.tsv").toString(), loci=root.resolve("loci.tsv").toString();
			final String graph=root.resolve("graph.tsv").toString(), exclude=root.resolve("exclude.txt").toString();
			final String fasta=root.resolve("contexts.fa").toString();
			NcrnaLearningPopulationTest.write(roster, NcrnaLearningPopulationTest.roster());
			NcrnaLearningPopulationTest.write(loci, NcrnaLearningPopulationTest.loci());
			NcrnaLearningPopulationTest.write(graph, NcrnaLearningPopulationTest.graph());
			NcrnaLearningPopulationTest.write(exclude, "gene0\ngene2\n");
			final NcrnaLearningPopulation population=NcrnaLearningPopulation.load(NcrnaSupportIO.read(graph),
				new String[]{"common", "rare"}, roster, loci, exclude, 42);
			checkContexts(population);
			NcrnaLearningPopulationTest.write(fasta, corpus(21, false));
			final boolean trim=Shared.TRIM_READ_DESCRIPTION, upper=Read.TO_UPPER_CASE;
			final NcrnaLearningTables train=NcrnaLearningTables.build(population, fasta, 3, 0, 2, 1, false, 3);
			final NcrnaLearningTables refit=NcrnaLearningTables.build(population, fasta, 3, 0, 2, 1, true, 2);
			assert(Shared.TRIM_READ_DESCRIPTION==trim && Read.TO_UPPER_CASE==upper)
				: "Context processing must restore its process-wide reader settings";
			assert(refit.records==21 && refit.learning==18 && train.records==21)
				: "All original contexts survive; only the eighteen nonpartial, nonexcluded loci may teach";
			assert(refit.start.counted()==15 && refit.start.missing()==2 && refit.start.ambiguous()==1
				&& refit.stop.counted()==16 && refit.stop.missing()==2 && refit.stop.ambiguous()==0)
				: "Refit has two missing spans at each end and one ambiguous start; missing GC still contributes to tables";
			assert(Arrays.equals(refit.availability, new long[]{18, 1, 1, 1}))
				: "Original population has eighteen available, one missing flank, one malformed span and one missing-GC context";
			assert(refit.statusByGene[1]==NcrnaLearningContext.NO_GC
				&& refit.statusByGene[5]==NcrnaLearningContext.NO_FLANK
				&& refit.statusByGene[6]==NcrnaLearningContext.MALFORMED_SPAN
				&& refit.statusByGene[20]==NcrnaLearningContext.AVAILABLE)
				: "Per-locus availability stays on the original axis, including partial loci, for vector-shard assignment";
			int teaching=0, missing=0, ambiguous=0;
			for(int gene=1; gene<20; gene++){
				if(gene==2 || population.split().validationByGene[gene]){continue;}
				teaching++;
				if(gene==5 || gene==6){missing++;}
				if(gene==3){ambiguous++;}
			}
			assert(train.learning==teaching && train.start.counted()==teaching-missing-ambiguous
				&& train.stop.counted()==teaching-missing)
				: "Tuning tables must use exactly the eligible TRAIN membership, excluding every validation locus";
			assert(refit.start.rawCount(bytes("AAA"), 0)==15 && refit.stop.rawCount(bytes("AAT"), 0)==16)
				: "The known context sequences produce literal AAA start and AAT stop counts";
			NcrnaLearningPopulationTest.write(fasta, corpus(20, false));
			fails(()->NcrnaLearningTables.build(population, fasta, 3, 0, 2, 1, true, 2));
			NcrnaLearningPopulationTest.write(fasta, corpus(21, true));
			fails(()->NcrnaLearningTables.build(population, fasta, 3, 0, 2, 1, true, 3));
			System.out.println("Learning contexts passed: complete streamed join, missing/duplicate rejection, metadata/geometry checks, training-only counts and all-eligible refit");
		}finally{
			Shared.setBufferLen(buffer);
		}
	}

	private static void checkContexts(final NcrnaLearningPopulation population){
		final NcrnaLearningContext c=new NcrnaLearningContext();
		final String header="gene1 source=genome1;lflank=20;rflank=30;contig_gc=0.5";
		c.set(read(header, sequence(1)), population);
		assert(c.gene==1 && c.start==20 && c.stop==119 && c.gc==0.5f && c.hasFeatures())
			: "Context coordinates are zero-based L through N-R-1, independent of header separators";
		c.set(read(header+" contig=c strand=+ start=21 stop=120 family=RF00174 status=FULL description", sequence(1)), population);
		assert(c.hasFeatures() && c.start==20 && c.stop==119 && c.gc==0.5f)
			: "Extractor provenance and unkeyed descriptions must not alter training geometry or GC";
		fails(()->c.set(read(header.replace("lflank=", "lflnk="), sequence(1)), population));
		fails(()->c.set(read(header+" mystery=1", sequence(1)), population));
		final ByteBuilder fixed=new ByteBuilder();
		for(final float expected:new float[]{0f, 1f, 0.2f, 0.33333334f, 0.0001f, 1f/Integer.MAX_VALUE}){
			final String value=fixed.clear().append((double)expected, 18).toString();
			assert(value.indexOf('e')<0 && value.indexOf('E')<0) : "Fixed-decimal GC must never use exponent notation: "+value;
			c.set(read(header.replace("contig_gc=0.5", "contig_gc="+value), sequence(1)), population);
			assert(Float.floatToIntBits(c.gc)==Float.floatToIntBits(expected))
				: "Serialized GC must preserve the caller's float feature, including a single GC base in an int-sized contig: "+value;
		}
		c.set(read(header.replace(";contig_gc=0.5", ""), sequence(1)), population);
		assert(c.hasGeometry() && !c.hasFeatures() && c.status==NcrnaLearningContext.NO_GC)
			: "A missing whole-contig GC prevents vectors but not endpoint counting";
		c.set(read(header.replace(";rflank=30", ""), sequence(1)), population);
		assert(!c.hasGeometry() && c.status==NcrnaLearningContext.NO_FLANK) : "No guessed flank geometry";
		c.set(read(header, new byte[10]), population);
		assert(c.status==NcrnaLearningContext.MALFORMED_SPAN) : "A short or impossible span is an explicit limitation";
		fails(()->c.set(read(header, new byte[151]), population));
		fails(()->c.set(read(header.replace("lflank=20", "lflank=21"), sequence(1)), population));
		fails(()->c.set(read(header.replace("lflank=20", "lflank=21").replace(";rflank=30", ""), sequence(1)), population));
		fails(()->c.set(read(header.replace("lflank=20", "lflank=2147483647"), sequence(1)), population));
		fails(()->c.set(read(header.replace("source=genome1", "source=genome0"), sequence(1)), population));
		fails(()->c.set(read(header+" source=genome1", sequence(1)), population));
		fails(()->c.set(read(header.replace("gene1 ", "unknown "), sequence(1)), population));
		fails(()->c.set(read(header.replace("source=", "xsource="), sequence(1)), population));
		fails(()->c.set(read(header.replace("lflank=20", "lflank=2147483648"), sequence(1)), population));
		fails(()->c.set(read(header.replace("contig_gc=0.5", "contig_gc=1.1"), sequence(1)), population));
		fails(()->c.set(read(header.replace("contig_gc=0.5", "contig_gc=garbage"), sequence(1)), population));
		fails(()->c.set(read(header.replace("contig_gc=0.5", "contig_gc=1.0E-4garbage"), sequence(1)), population));
		c.set(read(header.replace("contig_gc=0.5", "contig_gc=1.0E-4"), sequence(1)), population);
		assert(c.hasFeatures() && c.gc==0.0001f) : "Legacy exponent GC must be parsed as the complete fraction, not truncated to its mantissa";
		c.set(read(header.replace("contig_gc=0.5", "contig_gc=NaN"), sequence(1)), population);
		assert(c.status==NcrnaLearningContext.NO_GC) : "Explicit unavailable GC stays unavailable, never zero";
	}

	private static String corpus(final int count, final boolean duplicate){
		assert(count==20 || count==21) : "Fixture omission removes the last original, partial locus";
		final ByteBuilder b=new ByteBuilder();
		for(int gene=0; gene<count+(duplicate ? 1 : 0); gene++){
			final int id=gene==21 ? 0 : gene;
			b.append(">gene").append(id).append(" source=genome").append(id==20 ? 0 : id).append(" lflank=20");
			if(id!=5){b.append(" rflank=30");}
			if(id!=1){b.append(" contig_gc=0.5");}
			b.nl().append(sequence(id)).nl();
		}
		return b.toString();
	}

	private static byte[] sequence(final int gene){
		if(gene==6){return bytes("AAAAAAAAAA");}
		final byte[] bases=new byte[150];
		Arrays.fill(bases, 0, 20, (byte)'G');
		Arrays.fill(bases, 20, 120, (byte)'A');
		Arrays.fill(bases, 120, 150, (byte)'T');
		if(gene==3){bases[20]='N';}
		return bases;
	}

	private static Read read(final String header, final byte[] sequence){return new Read(sequence, null, header, 0, false);}
	private static byte[] bytes(final String value){return value.getBytes(StandardCharsets.US_ASCII);}

	private static void fails(final Runnable action){
		try{
			action.run();
		}catch(IllegalArgumentException expected){return;}
		throw new AssertionError("Invalid original context must fail instead of disappearing from the learning population");
	}
}
