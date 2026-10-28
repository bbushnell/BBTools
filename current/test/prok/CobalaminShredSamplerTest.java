package test.prok;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import map.ObjectMap;
import parse.LineParser1;
import prok.CobalaminShredSampler;
import prok.NcrnaCmLocusAudit;
import prok.NcrnaCorpusSample;
import prok.RrnaResourceIO;
import structures.ByteBuilder;

/** Independent extraction/mask checks, fixed-length conservation and source exclusion.
 * @author Raiden
 */
public final class CobalaminShredSamplerTest {

	public static void main(String[] args)throws Exception{
		check(args.length==1, "Supply a fresh fixture directory");final Path root=Paths.get(args[0]);Files.createDirectory(root);
		final String[] families={"RF00174", "RF01482", "RF01689"};final int[] starts={10, 70, 120}, stops={15, 80, 130};
		final byte[] bases=new byte[200];for(int i=0; i<bases.length; i++){bases[i]=(byte)"ACGT".charAt(i&3);}for(int i=44; i<50; i++){bases[i]='N';}
		final ByteBuilder roster=new ByteBuilder(NcrnaCorpusSample.GENOME_HEADER+'\n'), panel=new ByteBuilder();
		final ByteBuilder[] masks=new ByteBuilder[3];for(int f=0; f<3; f++){masks[f]=new ByteBuilder(NcrnaCmLocusAudit.HEADER+'\n');}
		for(int g=0; g<6; g++){
			final String name="g"+g, assembly=g==3?"GCF_0.2":"GCF_"+g+".1";final int taxid=g==2?1:g+1;
			final Path fasta=root.resolve(name+".fa.gz");write(fasta, ">contig\n"+new String(bases, StandardCharsets.US_ASCII)+'\n');
			roster.append(g).tab().append(name).tab().append(assembly).tab().append(taxid).append("\tbacteria\tP\t").append(root.resolve(name+".gff").toString())
				.tab().append(fasta.toString()).append("\t3\t200\t1\n");
			if(g<2){panel.append(g==0?"00000000000000000001":"00000000000000000002").tab().append(fasta.toString()).nl();}
			for(int f=0; f<3; f++){
				masks[f].append(families[f]).append('_').append(g).append("_0\t").append(name).tab().append(assembly).tab().append(taxid)
					.append("\tbacteria\tP\tcontig\t").append(f==1?"-":"+").tab().append(starts[f]).append('-').append(stops[f]).tab().append(stops[f]-starts[f]+1)
					.append("\t1\tx\tfalse\tfalse\tANNOTATION_FULL_COPY\t200\tfalse\t0\t100\t200\t1\n");
			}
		}
		write(root.resolve("genomes.tsv"), roster.toString());write(root.resolve("panel.tsv"), panel.toString());
		for(int f=0; f<3; f++){write(root.resolve(families[f]+".tsv"), masks[f].toString());}
		final ByteBuilder lengths=new ByteBuilder();for(int n:new int[]{17, 18, 20, 23}){lengths.append('>').append(n).nl();for(int i=0; i<n; i++){lengths.append('A');}lengths.nl();}
		write(root.resolve("lengths.fa"), lengths.toString());run(root, "panel", "a", "lengths.fa");run(root, "nonpanel", "b", "lengths.fa");run(root, "panel", "repeat", "lengths.fa");
		for(String suffix:new String[]{".shreds.tsv", ".genomes.tsv", ".lengths.tsv", ".summary.tsv"}){
			check(Arrays.equals(Files.readAllBytes(root.resolve("a"+suffix)), Files.readAllBytes(root.resolve("repeat"+suffix))), "Fixed seed and source order must reproduce selection: "+suffix);
		}
		check(Arrays.equals(Files.readAllBytes(root.resolve("a.lengths.tsv")), Files.readAllBytes(root.resolve("b.lengths.tsv"))), "Both arms must use the identical empirical-length schedule");
		inspect(root, "a", "panel", bases, starts, stops);inspect(root, "b", "nonpanel", bases, starts, stops);
		final ByteBuilder impossible=new ByteBuilder();for(int n=0; n<4; n++){impossible.append('>').append(n).nl();for(int i=0; i<500; i++){impossible.append('A');}impossible.nl();}
		write(root.resolve("impossible.fa"), impossible.toString());boolean failed=false;
		try{run(root, "panel", "impossible", "impossible.fa");}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "No eligible start must fail instead of silently redrawing lengths or shrinking the denominator");
		write(root.resolve("RF01482.tsv"), masks[1].toString().replace("\tcontig\t", "\tmissing\t"));failed=false;
		try{run(root, "panel", "missing", "lengths.fa");}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "A mask absent from source FASTA must not silently become unmasked background");
		System.out.println("COBALAMIN_SHRED_SAMPLER_TEST_PASS");
	}
	private static void inspect(Path root, String prefix, String arm, byte[] bases, int[] starts, int[] stops){
		final ObjectMap<String,byte[]> cores=new ObjectMap<String,byte[]>(String.class, byte[].class), contexts=new ObjectMap<String,byte[]>(String.class, byte[].class);
		RrnaResourceIO.stream(root.resolve(prefix+".cores.fa.gz").toString(), r->{check(cores.put(r.id, r.bases)==null, "Core IDs must be unique");});
		RrnaResourceIO.stream(root.resolve(prefix+".contexts.fa.gz").toString(), r->{check(contexts.put(r.id, r.bases)==null, "Context IDs must be unique");});
		final ByteFile in=RrnaResourceIO.open(root.resolve(prefix+".shreds.tsv").toString());final LineParser1 p=new LineParser1('\t');int count=0, clipped=0;
		try{RrnaResourceIO.header(in, CobalaminShredSampler.HEADER, "shreds");for(byte[] row; (row=in.nextLine())!=null; ){
			p.set(row);count++;check(p.terms()==13 && p.termEquals(arm, 1), "Every selected core needs exact arm and coordinate provenance");
			final String id=p.parseString(0), genome=p.parseString(2);check(arm.equals("panel")?genome.equals("g0") || genome.equals("g1"):genome.equals("g4") || genome.equals("g5"),
				"Nonpanel must exclude panel paths, shared taxids and other versions of panel assemblies");
			final int lo=p.parseInt(6), hi=p.parseInt(7), left=p.parseInt(8), right=p.parseInt(9), length=p.parseInt(12);
			check(length==hi-lo+1 && (length==17 || length==18 || length==20 || length==23), "Core length must remain a draw from the empirical positive distribution");
			check(left==Math.max(1, lo-30) && right==Math.min(200, hi+30) && p.parseInt(10)==lo-left+1 && p.parseInt(11)==hi-left+1,
				"Padded contexts must clip at real contig edges and preserve core coordinates for CM-overlap exclusion");
			if(left==1 || right==200){clipped++;}
			for(int i=0; i<3; i++){check(hi<starts[i] || lo>stops[i], "Both strands and all three source families must mask background");}
			check(hi<45 || lo>50, "Ambiguous bases must not enter core background");
			check(Arrays.equals(cores.get(id), Arrays.copyOfRange(bases, lo-1, hi)) && Arrays.equals(contexts.get(id), Arrays.copyOfRange(bases, left-1, right)),
				"Native FASTA output must match independent source-coordinate extraction");
		}}finally{check(!in.close(), "Shred manifest read failed");}
		check(count==200 && cores.size()==200 && contexts.size()==200 && clipped>0, "All fixed schedule entries and clipped-edge cases must survive");
	}
	private static void run(Path root, String arm, String out, String lengths){
		CobalaminShredSampler.main(new String[]{"genomes="+root.resolve("genomes.tsv"), "panel="+root.resolve("panel.tsv"), "lengths="+root.resolve(lengths),
			"loci="+root.resolve("RF00174.tsv")+","+root.resolve("RF01482.tsv")+","+root.resolve("RF01689.tsv"), "outprefix="+root.resolve(out), "arm="+arm,
			"count=200", "genomecount=2", "expectedgenomes=6", "expectedpanel=2", "expectedlengths=4", "expectedloci=18", "seed=20261026", "flank=30"});
	}
	/** Negative cases deliberately replace fixture inputs inside the fresh test root. */
	private static void write(Path path, String text){
		final ByteStreamWriter out=new ByteStreamWriter(FileFormat.testOutput(path.toString(), FileFormat.TEXT, null, true, true, false, false));out.start();
		try{out.print(text);}finally{RrnaResourceIO.finish(out);}
	}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
