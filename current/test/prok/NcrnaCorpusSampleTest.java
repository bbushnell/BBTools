package test.prok;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Locale;

import fileIO.ByteFile;
import parse.LineParser1;
import prok.NcrnaCmLocusAudit;
import prok.NcrnaCorpusSample;

/** Transitive source exclusion and row-order-independent sampling fixtures.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaCorpusSampleTest {
	public static void main(String[] args)throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Supply an unused fixture directory");}
		final Path root=Paths.get(args[0]);Files.createDirectory(root);
		final ArrayList<String> genomes=new ArrayList<String>(), loci=new ArrayList<String>(), hashes=new ArrayList<String>();
		for(int i=0;i<12;i++){
			final String name=String.format(Locale.ROOT, "g%02d", i), assembly=i==2?"GCF_2.2":"GCF_"+(i+1)+".1";
			final int taxid=i==1?1:i+1, length=i==0?100:i==8?80:i==9?450:190, start=i==10?1:1001, left=i==10?0:200;
			final String phylum=i<7?"P":"Q", fasta="/fixture/"+name+".fna.gz";
			genomes.add(i+"\t"+name+'\t'+assembly+'\t'+taxid+"\tbacteria\t"+phylum+"\t/fixture/"+name+".gff.gz\t"+fasta+"\t1\t2000\t1\n");
			hashes.add(String.format(Locale.ROOT, "%020x", i==3?3:i+1)+'\t'+fasta+'\n');
			loci.add("locus_"+name+'\t'+name+'\t'+assembly+'\t'+taxid+"\tbacteria\t"+phylum+"\tcontig\t+\t"+start+'-'+(start+length-1)+'\t'+length
				+"\t1\tfeature\t"+(i==7)+"\tfalse\t"+(i==7?"SOURCE_PARTIAL":"ANNOTATION_FULL_COPY")+"\t2000\t"+(i==10)+'\t'+(i==11?1:0)+'\t'+left+"\t200\t1\n");
		}
		write(root.resolve("genomes.tsv"), NcrnaCorpusSample.GENOME_HEADER+'\n', genomes);
		write(root.resolve("loci.tsv"), NcrnaCmLocusAudit.HEADER+'\n', loci);
		write(root.resolve("hashes.tsv"), "", hashes);
		write(root.resolve("pilot.tsv"), NcrnaCorpusSample.GENOME_HEADER+'\n', new ArrayList<String>(Collections.singleton(genomes.get(0))));
		run(root, "a", 2);
		final java.util.List<String> eligible=Files.readAllLines(root.resolve("a.eligible.ids.txt"), StandardCharsets.US_ASCII);
		check(eligible.size()==7 && eligible.contains("locus_g00") && eligible.contains("locus_g03") && eligible.contains("locus_g06"),
			"Full-population export must include eligible pilot-linked groups and cannot be replaced by either validation sample arm");
		final ByteFile in=ByteFile.makeByteFile(root.resolve("a.loci.tsv").toString(), false);final LineParser1 p=new LineParser1('\t');
		int probability=0, challenge=0, edgePriority=0;
		try{
			check(in.nextLine()!=null, "Sample ledger header required");
			for(byte[] row;(row=in.nextLine())!=null;){
				p.set(row);check(p.terms()==29, "Every selected row keeps source metadata and sampling provenance");
				check(p.termEquals("NOT_ASSIGNED", 26), "Sampling must not invent a train/holdout split");
				if(p.termEquals("PROBABILITY", 25)){
					probability++;check(p.termEquals("false", 22) && p.termEquals("true", 23), "Probability arm cannot contain pilot-linked or ineligible sources");
					final String genome=p.parseString(1);check(genome.equals("g04") || genome.equals("g05") || genome.equals("g06"), "Transitive taxid/accession/hash group or flags leaked into probability arm: "+genome);
				}else{
					check(p.termEquals("CHALLENGE", 25), "Unknown sample arm");challenge++;
					if(p.termEquals("P", 5)){check(p.termEquals("g00", 1) && p.parseInt(28)==0, "Band-edge challenge priority must beat interior loci within the same stratum");edgePriority++;}
				}
			}
		}finally{check(!in.close(), "Sample output read failed");}
		check(probability==2 && challenge==6, "Sample arms must retain requested distinct counts");
		check(edgePriority==1, "Fixture must test a boundary-prioritized challenge choice");
		final ByteFile groups=ByteFile.makeByteFile(root.resolve("a.groups.tsv").toString(), false);int linked=0;
		try{groups.nextLine();for(byte[] row;(row=groups.nextLine())!=null;){
			p.set(row);if(p.termEquals("true", 5)){check(p.termEquals("g00", 4), "Transitive pilot component must have canonical group identity");linked++;}
		}}finally{check(!groups.close(), "Group table read failed");}
		check(linked==4, "Pilot exclusion must propagate through taxid, assembly version, and identical FASTA hash");
		final ByteFile strata=ByteFile.makeByteFile(root.resolve("a.strata.tsv").toString(), false);int excluded=0;
		try{strata.nextLine();for(byte[] row;(row=strata.nextLine())!=null;){
			p.set(row);check(p.terms()==9 && p.termEquals("PROXY_ELIGIBLE_AND_NONPILOT_GROUP", 8), "Sampling fraction must name its conditional eligible population");
			check(p.parseInt(7)==3, "Probability denominator is the untouched eligible pool, not total stratum population");
			if(p.parseInt(4)>0){check(p.parseInt(5)>0, "Every nonempty challenge stratum must receive a sample");}
			if(p.parseInt(2)==0){check(p.parseInt(6)==0 && p.parseInt(3)==0, "Excluded strata must not claim nonzero probability inclusion");excluded++;}
			else{check(p.parseInt(6)==2, "Uniform eligible-locus probability is requested sample divided by eligible population");}
		}}finally{check(!strata.close(), "Stratum table read failed");}
		check(excluded>0, "Fixture must exercise zero-availability weighting");
		Collections.reverse(genomes);Collections.reverse(loci);Collections.reverse(hashes);
		write(root.resolve("genomes.tsv"), NcrnaCorpusSample.GENOME_HEADER+'\n', genomes);
		write(root.resolve("loci.tsv"), NcrnaCmLocusAudit.HEADER+'\n', loci);write(root.resolve("hashes.tsv"), "", hashes);
		run(root, "b", 2);
		for(String suffix:new String[]{".ids.txt", ".eligible.ids.txt", ".loci.tsv", ".strata.tsv", ".groups.tsv", ".summary.tsv"}){
			check(Arrays.equals(Files.readAllBytes(root.resolve("a"+suffix)), Files.readAllBytes(root.resolve("b"+suffix))), "Input row order changed seeded sample: "+suffix);
		}
		boolean failed=false;try{run(root, "insufficient", 4);}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "Insufficient untouched probability pool must fail rather than import pilot groups");
		NcrnaCorpusSample.main(new String[]{"loci="+root.resolve("loci.tsv"), "genomes="+root.resolve("genomes.tsv"), "hashes="+root.resolve("hashes.tsv"),
			"pilotgenomes="+root.resolve("pilot.tsv"), "outprefix="+root.resolve("census"), "expected=12", "expectedgenomes=12", "minlen=500", "maxlen=600",
			"probability=0", "challenge=0", "seed=20261025"});
		final ByteFile census=ByteFile.makeByteFile(root.resolve("census.loci.tsv").toString(), false);int censusRows=0;
		try{
			check(census.nextLine()!=null, "Census preserves the metadata schema even when the eligible population is empty");
			for(byte[] row; (row=census.nextLine())!=null; ){
				p.set(row);censusRows++;
				check(p.terms()==29 && p.termEquals("false", 23) && p.termEquals("CENSUS", 25) && p.termEquals("NOT_ASSIGNED", 26),
					"Census must retain every excluded source without pretending to sample or split it");
			}
		}finally{check(!census.close(), "Census output read failed");}
		check(censusRows==12 && Files.size(root.resolve("census.eligible.ids.txt"))==0, "Zero eligible members must not erase the twelve-row source census");
		System.out.println("NCRNA_CORPUS_SAMPLE_TEST_PASS");
	}
	private static void run(Path root, String name, int probability){
		NcrnaCorpusSample.main(new String[]{"loci="+root.resolve("loci.tsv"), "genomes="+root.resolve("genomes.tsv"), "hashes="+root.resolve("hashes.tsv"),
			"pilotgenomes="+root.resolve("pilot.tsv"), "outprefix="+root.resolve(name), "expected=12", "expectedgenomes=12", "minlen=100", "maxlen=300",
			"probability="+probability, "challenge=6", "seed=20261025"});
	}
	private static void write(Path path, String header, ArrayList<String> rows)throws Exception{
		final StringBuilder b=new StringBuilder(header);for(String row:rows){b.append(row);}Files.write(path, b.toString().getBytes(StandardCharsets.US_ASCII));
	}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
