package test.prok;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.LineParser1;
import prok.NcrnaClanSeedPrep;
import prok.NcrnaCmLocusAudit;
import prok.NcrnaCorpusSample;
import prok.RrnaResourceIO;
import structures.ByteBuilder;

/** Checks sibling-only source hashes, immutable raw rows and corrected alias geometry.
 * @author Brian, Raiden
 */
public final class NcrnaClanSeedPrepTest {

	public static void main(String[] args)throws Exception{
		check(args.length==1, "Supply a fresh fixture directory");final Path root=Paths.get(args[0]);Files.createDirectory(root);
		final String[] families={"RF00174", "RF01482", "RF01689"};
		final int[][] present={{1, 1, 0}, {1, 0, 1}, {0, 1, 1}};
		final ByteBuilder manifest=new ByteBuilder("family\tloci\tgenomes\thashes\texpected\n");
		final ArrayList<String> raw=new ArrayList<String>();String secondRoster=null, secondHashes=null;
		for(int f=0; f<3; f++){
			final ByteBuilder genomes=new ByteBuilder(NcrnaCorpusSample.GENOME_HEADER+'\n'), loci=new ByteBuilder(NcrnaCmLocusAudit.HEADER+'\n'), hashes=new ByteBuilder();
			for(int g=0; g<3; g++){
				final String name="g"+g, phylum=g<2?"A":"B";final boolean has=present[f][g]>0;
				genomes.append(g).tab().append(name).append("\tGCF_").append(g).append(".1\t").append(g+1).append("\tbacteria\t").append(phylum)
					.tab().append(root.resolve(name+".gff.gz").toString()).tab().append(root.resolve(name+".fna.gz").toString())
					.tab().append(has?1:0).tab().append(has?1000:0).tab().append(has?1:0).nl();
				if(has){
					hashes.append(String.format(Locale.ROOT, "%020x", g+1)).tab().append(root.resolve(name+".fna.gz").toString()).nl();
					final int start=g==1?201:101, length=f==0 && g==0?20:30;
					final String row=families[f]+"_"+g+"_0\t"+name+"\tGCF_"+g+".1\t"+(g+1)+"\tbacteria\t"+phylum+
						"\tcontig\t+\t"+start+"-"+(start+length-1)+"\t"+length+"\t1\tx\tfalse\tfalse\tANNOTATION_FULL_COPY\t1000\tfalse\t0\t100\t100\t1";
					loci.append(row).nl();raw.add(row);
				}
			}
			write(root.resolve(families[f]+".genomes.tsv"), genomes.toString());write(root.resolve(families[f]+".loci.tsv"), loci.toString());
			write(root.resolve(families[f]+".hashes.tsv"), hashes.toString());
			if(f==1){secondRoster=genomes.toString();secondHashes=hashes.toString();}
			manifest.append(families[f]).tab().append(root.resolve(families[f]+".loci.tsv").toString()).tab().append(root.resolve(families[f]+".genomes.tsv").toString())
				.tab().append(root.resolve(families[f]+".hashes.tsv").toString()).append("\t2\n");
		}
		write(root.resolve("cohorts.tsv"), manifest.toString());
		final String input=">RF00174_0_0 source=g0 member_family=RF00174 build_group=RF00174\n"+bases('A', 20)+
			"\n>RF00174_1_0 source=g1 contig=contig strand=+ start=201 stop=219 member_family=RF00174 build_group=RF00174\n"+bases('G', 19)+
			"\n>RF01482_2_0 source=g2 source_family=RF01482 contig=contig strand=+ start=99 stop=125 member_family=RF01482 build_group=RF01482\n"+bases('C', 27)+
			"\n>RF01689_1_0 source=g1 source_family=RF01689 contig=contig strand=+ start=201 stop=219 member_family=B12clade build_group=RF00174\n"+bases('G', 19)+'\n';
		write(root.resolve("input.fa"), input);run(root, "input.fa", "good");
		final ByteFile original=RrnaResourceIO.open(root.resolve("good.raw.loci.tsv.gz").toString());
		try{RrnaResourceIO.header(original, NcrnaCmLocusAudit.HEADER, "raw");for(String row:raw){check(Arrays.equals(original.nextLine(), row.getBytes(StandardCharsets.UTF_8)), "Original annotation rows must remain byte-identical and complete");}check(original.nextLine()==null, "No invented raw source rows");}
		finally{check(!original.close(), "Raw fixture read failed");}
		final LineParser1 p=new LineParser1('\t');final ByteFile roster=RrnaResourceIO.open(root.resolve("good.genomes.tsv").toString());
		try{roster.nextLine();for(int g=0; g<3; g++){p.set(roster.nextLine());check(p.parseInt(8)==2 && p.parseLong(9)==1000 && p.parseInt(10)==1,
			"Merged counts sum all families while sibling-only positives supply actual FASTA size observations");}check(roster.nextLine()==null, "Genome roster must not gain rows");}
		finally{check(!roster.close(), "Genome fixture read failed");}
		final ByteFile hashes=RrnaResourceIO.open(root.resolve("good.hashes.sha80").toString());int hashCount=0;
		try{for(byte[] row; (row=hashes.nextLine())!=null; ){p.set(row);hashCount++;check(p.terms()==2, "Hash provenance stays sha80/path");}}
		finally{check(!hashes.close(), "Hash fixture read failed");}check(hashCount==3, "Sibling-only positive genome needs its hash; shared paths deduplicate");
		final ByteFile seed=RrnaResourceIO.open(root.resolve("good.seed.loci.tsv").toString());
		try{
			RrnaResourceIO.header(seed, NcrnaClanSeedPrep.SEED_HEADER, "seed");p.set(seed.nextLine());check(p.terms()==12 && p.termEquals("101-120", 8), "Uncorrected proxy retains its original interval");
			p.set(seed.nextLine());check(p.termEquals("201-219", 8) && p.termEquals("RF00174", 10), "Recovered core uses corrected geometry");
			p.set(seed.nextLine());check(p.termEquals("99-125", 8) && p.termEquals("RF01482", 9), "Corrected core may extend beyond its original annotation");
			p.set(seed.nextLine());check(p.termEquals("201-219", 8) && p.termEquals("RF01689", 9) && p.termEquals("B12clade", 10)
				&& p.termEquals("RF00174", 11) && seed.nextLine()==null, "Coincident cores must preserve distinct source, biological and operational labels");
		}finally{check(!seed.close(), "Seed fixture read failed");}
		write(root.resolve("RF01482.hashes.tsv"), secondHashes.replace("00000000000000000001", "00000000000000000009"));reject(root, "input.fa", "conflicting_hash");
		write(root.resolve("RF01482.hashes.tsv"), secondHashes);
		write(root.resolve("RF01482.genomes.tsv"), secondRoster.replace("\tbacteria\tB\t", "\tbacteria\tCHANGED\t"));reject(root, "input.fa", "taxonomy");
		write(root.resolve("RF01482.genomes.tsv"), secondRoster);
		write(root.resolve("bad_bounds.fa"), input.replace("start=99 stop=125", "start=99 stop=126"));reject(root, "bad_bounds.fa", "bounds");
		write(root.resolve("missing.fa"), input.substring(0, input.indexOf(">RF01689")));reject(root, "missing.fa", "missing");
		write(root.resolve("unknown.fa"), input.replace(">RF00174_0_0 ", ">RF00174_999_0 "));reject(root, "unknown.fa", "unknown");
		System.out.println("NCRNA_CLAN_SEED_PREP_TEST_PASS");
	}
	private static void run(Path root, String input, String out){
		NcrnaClanSeedPrep.main(new String[]{"in="+root.resolve(input), "cohorts="+root.resolve("cohorts.tsv"), "outprefix="+root.resolve(out), "expected=4", "expectedsource=6", "expectedgenomes=3"});
	}
	private static void reject(Path root, String input, String out){
		boolean failed=false;try{run(root, input, out);}catch(IllegalArgumentException expected){failed=true;}
		check(failed && !Files.exists(root.resolve(out+".raw.loci.tsv.gz")), "Invalid source or admitted membership must fail before publishing: "+out);
	}
	private static String bases(char base, int count){final ByteBuilder b=new ByteBuilder();for(int i=0; i<count; i++){b.append(base);}return b.toString();}
	/** Negative cases deliberately replace fixture inputs inside the fresh test root. */
	private static void write(Path path, String text){
		final ByteStreamWriter out=new ByteStreamWriter(FileFormat.testOutput(path.toString(), FileFormat.TEXT, null, true, true, false, false));out.start();
		try{out.print(text);}finally{RrnaResourceIO.finish(out);}
	}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
