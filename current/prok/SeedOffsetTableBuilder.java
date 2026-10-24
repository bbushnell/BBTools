package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import map.LongHashSet;
import map.LongObjectMap;
import parse.Parser;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.ReadInputStream;
import structures.ByteBuilder;

/** Builds a seed endpoint-offset table from annotation-bounded, sense-oriented genes.
 * The orchestrator supplies the sha80 of the canonical sorted caller-visible seed set;
 * the loader later checks that value against the paired seed resource.
 */
public final class SeedOffsetTableBuilder {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, SeedOffsetTableBuilder.class, false);
		args=pp.args;
		final Parser parser=new Parser();
		parser.overwrite=false;
		String in=null, seeds=null, out=null, seedSha80=null, trainingSource="genomecut";
		int k=17, expectedGenes=-1;
		for(String arg : args){
			final int eq=arg.indexOf('=');
			final String key=(eq<0 ? arg : arg.substring(0, eq)).toLowerCase(Locale.ROOT);
			final String value=(eq<0 ? null : arg.substring(eq+1));
			if(key.equals("seeds")){seeds=value;}
			else if(key.equals("seedsha80") || key.equals("seedsetsha80")){seedSha80=value;}
			else if(key.equals("trainingsource")){trainingSource=value.toLowerCase(Locale.ROOT);}
			else if(key.equals("expectedgenes")){expectedGenes=Integer.parseInt(value);}
			else if(key.equals("k")){k=Integer.parseInt(value);}
			else if(parser.parse(arg, key, value)){}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		in=parser.in1;
		out=parser.out1;
		if(in==null || seeds==null || out==null || seedSha80==null){
			throw new IllegalArgumentException("Required: in= seeds= out= seedsha80= [k=17] [trainingsource=genomecut|genome_id|annotated] [expectedgenes=] [overwrite=f]");
		}
		if(expectedGenes<1){throw new IllegalArgumentException("Real seed-offset builds require expectedgenes>0");}
		if(trainingSource.equals("all")){throw new IllegalArgumentException("trainingsource=all is test-only; use genomecut or genome_id");}
		build(in, seeds, out, seedSha80, k, parser.overwrite, trainingSource, expectedGenes, parser.maxReads);
	}

	public static void build(String in, String seedPath, String out, String seedSha80, int k, boolean overwrite){
		build(in, seedPath, out, seedSha80, k, overwrite, "all", -1);
	}

	public static void build(String in, String seedPath, String out, String seedSha80, int k, boolean overwrite,
			String trainingSource, int expectedGenes){
		build(in, seedPath, out, seedSha80, k, overwrite, trainingSource, expectedGenes, -1);
	}

	static void build(String in, String seedPath, String out, String seedSha80, int k, boolean overwrite,
			String trainingSource, int expectedGenes, long maxReads){
		if(k!=17){throw new IllegalArgumentException("Seed-offset resources require K=17, found "+k);}
		if(!trainingSource.equals("genomecut") && !trainingSource.equals("genome_id") && !trainingSource.equals("annotated") && !trainingSource.equals("all")){
			throw new IllegalArgumentException("Unknown trainingSource: "+trainingSource);
		}
		if(!SeedOffsetTable.validSha80(seedSha80)){throw new IllegalArgumentException("Invalid seed-set sha80: "+seedSha80);}
		final LongHashSet seeds=ProkObject.loadLongKmers(seedPath, k);
		if(seeds==null || seeds.size()<1){throw new IllegalArgumentException("Empty caller-visible seed set: "+seedPath);}
		final long[] sortedSeeds=seeds.toArray();
		Arrays.sort(sortedSeeds);
		final LongObjectMap<Moments> stats=new LongObjectMap<Moments>(Math.max(16, seeds.size()*2), Moments.class);
		for(long seed : sortedSeeds){stats.put(seed, new Moments());}

		if(trainingSource.equals("annotated")){Shared.TRIM_READ_DESCRIPTION=false;}
		final ArrayList<Read> genes=ReadInputStream.toReads(in, FileFormat.FA, maxReads);
		if(genes==null || genes.isEmpty()){throw new IllegalArgumentException("No training genes: "+in);}
		int retained=0;
		for(Read gene : genes){
			if(!retain(gene.id, trainingSource)){continue;}
			if(trainingSource.equals("genomecut")){validateGenomeCut(gene);}
			if(trainingSource.equals("annotated")){validateAnnotated(gene);}
			addGene(gene.bases, k, seeds, stats); retained++;
		}
		if(retained<1 || (expectedGenes>=0 && retained!=expectedGenes)){
			throw new IllegalStateException("Training-gene count mismatch: observed="+retained+", expected="+expectedGenes);
		}

		final ByteStreamWriter bsw=new ByteStreamWriter(FileFormat.testOutput(out, FileFormat.TEXT, null, true, overwrite, false, false));
		final ByteBuilder bb=new ByteBuilder();
		bsw.start();
		bsw.println("#SeedOffsetTable\t"+SeedOffsetTable.VERSION);
		bsw.println("#K\t"+k);
		bsw.println("#SeedCount\t"+sortedSeeds.length);
		bsw.println("#TrainingGenes\t"+retained);
		bsw.println("#SeedSetSha80\t"+seedSha80);
		bsw.println(SeedOffsetTable.HEADER);
		for(long seed : sortedSeeds){
			final Moments m=stats.get(seed);
			final long left10=(m.count==0 ? 0 : Math.round(m.leftMean*10));
			final long right10=(m.count==0 ? 0 : Math.round(m.rightMean*10));
			final long leftSd10=(m.count==0 ? 0 : Math.round(Math.sqrt(Math.max(0, m.leftM2/m.count))*10));
			final long rightSd10=(m.count==0 ? 0 : Math.round(Math.sqrt(Math.max(0, m.rightM2/m.count))*10));
			bb.clear().appendA48(seed).tab().appendA48(left10).tab().appendA48(right10).tab()
				.appendA48(leftSd10).tab().appendA48(rightSd10).tab().appendA48(m.count).nl();
			bsw.print(bb.toBytes());
		}
		if(bsw.poisonAndWait()){throw new RuntimeException("Write error: "+out);}
	}

	private static boolean retain(String header, String source){
		if(source.equals("all") || source.equals("annotated")){return true;}
		if(source.equals("genomecut")){
			final int space=header.indexOf(' ');
			return !(space>=0 && space+1<header.length() &&
				(header.startsWith("Bacteria;",space+1) || header.startsWith("Archaea;",space+1)
					|| header.startsWith("Eukaryota;",space+1)));
		}
		final boolean genome=header.indexOf("|ID=")>=0, silva=header.indexOf("|lcl|LSU_")>=0;
		if(genome==silva){throw new IllegalArgumentException("Expected exactly one 23S source marker: "+header);}
		return genome;
	}

	private static void validateGenomeCut(Read gene){
		final String h=gene.id; final int last=h.lastIndexOf('\t'), prev=(last<0 ? -1 : h.lastIndexOf('\t',last-1));
		if(prev<0 || last!=prev+2){throw new IllegalArgumentException("Genome-cut header lacks strand/coordinates: "+h);}
		final char strand=h.charAt(prev+1); if(strand!='+' && strand!='-'){throw new IllegalArgumentException("Invalid strand: "+h);}
		final String coord=h.substring(last+1); final int dash=coord.indexOf('-');
		if(dash<1){throw new IllegalArgumentException("Invalid coordinates: "+h);}
		final long start=Long.parseLong(coord.substring(0,dash)), stop=Long.parseLong(coord.substring(dash+1));
		if(start<0 || stop<start || stop-start+1!=gene.length()){
			throw new IllegalArgumentException("Genome-cut span/sequence mismatch: "+h+" length="+gene.length());
		}
	}

	/** CM-extracted genes retain their original 1-based inclusive genomic metadata;
	 * bases are already sense-oriented by the extraction producer. */
	private static void validateAnnotated(Read gene){
		final String h=gene.id;final int startIndex=fieldStart(h," start="),stopIndex=fieldStart(h," stop="),strandIndex=fieldStart(h," strand=");
		final long start=Parse.parseLong(h,startIndex,fieldEnd(h,startIndex)),stop=Parse.parseLong(h,stopIndex,fieldEnd(h,stopIndex));
		final char strand=h.charAt(strandIndex);
		if(start<1 || stop<start || stop-start+1!=gene.length() || fieldEnd(h,strandIndex)!=strandIndex+1 || (strand!='+' && strand!='-')){
			throw new IllegalArgumentException("Annotated gene span/strand does not match exact sense-oriented bases: "+h);
		}
		fieldStart(h," seqid=");
		fieldStart(h,h.contains(" assembly=") ? " assembly=" : " genome=");
	}
	private static int fieldStart(String header,String token){
		final int at=header.indexOf(token);
		if(at<0 || header.indexOf(token,at+1)>=0){throw new IllegalArgumentException("Missing or duplicate annotated field"+token+" in "+header);}
		final int start=at+token.length();
		if(start>=header.length() || Character.isWhitespace(header.charAt(start))){throw new IllegalArgumentException("Empty annotated field"+token+" in "+header);}
		return start;
	}
	private static int fieldEnd(String header,int start){
		assert(start>=0 && start<header.length()) : "Validated annotation field must begin inside its header";
		int end=start;while(end<header.length()&&!Character.isWhitespace(header.charAt(end))){end++;}return end;
	}

	static void addGene(byte[] bases, int k, LongHashSet seeds, LongObjectMap<Moments> stats){
		if(bases==null || bases.length<k){return;}
		final byte[] bton=AminoAcid.baseToNumber;
		final long mask=~((-1L)<<(2*k));
		long kmer=0;
		int len=0;
		for(int end=0; end<bases.length; end++){
			final int x=bton[bases[end]];
			if(x<0){len=0; kmer=0; continue;}
			kmer=((kmer<<2)|x)&mask;
			len++;
			if(len<k || !seeds.contains(kmer)){continue;}
			final int start=end-k+1;
			stats.get(kmer).add(start, bases.length-end-1);
		}
	}

	static final class Moments {
		void add(double left, double right){
			count++;
			double delta=left-leftMean; leftMean+=delta/count; leftM2+=delta*(left-leftMean);
			delta=right-rightMean; rightMean+=delta/count; rightM2+=delta*(right-rightMean);
		}
		long count;
		double leftMean, rightMean, leftM2, rightM2;
	}
}
