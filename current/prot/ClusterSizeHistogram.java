package prot;

import java.io.File;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.Parser;
import parse.PreParser;
import shared.Shared;

/**
 * Counts a representative-grouped MMseqs membership table and reports exact
 * membership-weighted coverage. The representative FASTA provides an independent
 * family inventory: every ID must occur in exactly one contiguous table group,
 * and every group must contain its representative exactly once as a member.
 * Member IDs are not globally deduplicated; each validated TSV row is one member.
 * Memory holds representative IDs and a primitive size-frequency array, not all
 * protein rows. Files are written only after input and count validation succeeds.
 * @author Brian Bushnell, Keqing
 */
public final class ClusterSizeHistogram {

	public static void main(String[] args){
		ClusterSizeHistogram x=new ClusterSizeHistogram(args);
		try{x.process();}finally{Shared.closeStream(x.log);}
	}

	private ClusterSizeHistogram(String[] args){
		PreParser pp=new PreParser(args, getClass(), false);
		args=pp.args;
		log=pp.outstream;
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			final int eq=args[i].indexOf('=');
			final String key=(eq<0 ? args[i] : args[i].substring(0, eq)).toLowerCase(Locale.ROOT);
			final String value=(eq<0 ? null : args[i].substring(eq+1));
			if(key.equals("reps")){reps=value;}
			else if(key.equals("expectedfamilies")){expectedFamilies=Long.parseLong(value);}
			else if(key.equals("expectedgenes")){expectedGenes=Long.parseLong(value);}
			else if(!parser.parse(args[i], key, value)){throw new IllegalArgumentException("Unknown argument: "+args[i]);}
		}
		in=parser.in1;
		out=parser.out1;
		require(in!=null && reps!=null && out!=null, "Required: in=membership.tsv reps=representatives.faa out=prefix");
		require(!parser.append, "Append would combine different census populations");
		require(expectedFamilies>=-1 && expectedGenes>=-1, "Expected counts must be nonnegative or -1 (unspecified)");
		require(new File(in).isFile() && new File(reps).isFile(), "Both inputs must be existing files");
		for(String suffix : SUFFIXES){require(!new File(out+suffix).exists(), "Refusing existing output: "+out+suffix);}
		Shared.setThreads(1);
		ByteFile.FORCE_MODE_BF1=true;
	}

	/** Removes IDs from the independent inventory as groups are encountered. */
	private void process(){
		final HashSet<String> unseen=readRepresentatives();
		final long inventory=unseen.size();
		ByteFile bf=ByteFile.makeByteFile(in, true);
		byte[] previous=null;
		int size=0, self=0;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				rows++;
				int tab=-1;
				for(int i=0; i<line.length; i++){
					if(line[i]=='\t'){
						if(tab>=0){throw new IllegalArgumentException("Extra TSV field at row "+rows);}
						tab=i;
					}else if(line[i]<=32){throw new IllegalArgumentException("Whitespace in ID at row "+rows);}
				}
				if(tab<=0 || tab>=line.length-1){throw new IllegalArgumentException("Expected rep<TAB>member at row "+rows);}
				if(!same(line, 0, tab, previous)){
					if(previous!=null){finishGroup(size, self);}
					previous=Arrays.copyOf(line, tab);
					final String id=new String(previous, StandardCharsets.US_ASCII);
					if(!unseen.remove(id)){throw new IllegalArgumentException("Unknown or noncontiguous repeated representative at row "+rows+": "+id);}
					size=0;
					self=0;
				}
				if(size>=Integer.MAX_VALUE-1){throw new IllegalArgumentException("Cluster size exceeds supported int histogram index at row "+rows);}
				size++;
				if(same(line, tab+1, line.length-tab-1, previous)){self++;}
				if(rows%10000000==0){log.println("Rows processed: "+rows+"; completed families: "+families);}
			}
		}finally{require(!bf.close(), "Membership reader reported an I/O failure");}
		if(previous!=null){finishGroup(size, self);}
		require(rows>0 && unseen.isEmpty() && families==inventory,
			"Incomplete census: rows="+rows+", groups="+families+", representatives="+inventory+", unseen="+unseen.size());
		require(expectedGenes<0 || rows==expectedGenes, "Member count differs: observed="+rows+", expected="+expectedGenes);
		require(expectedFamilies<0 || families==expectedFamilies,
			"Family count differs: observed="+families+", expected="+expectedFamilies);
		writeReports();
		log.println("CLUSTER_CENSUS_PASS families="+families+" proteins="+rows+" largest="+maxSize+" singletons="+hist[1]);
	}

	/** Reads FASTA header tokens, rejects duplicate IDs and empty records. */
	private HashSet<String> readRepresentatives(){
		final int capacity=(int)Math.min(1<<29, expectedFamilies>0 ? expectedFamilies*4/3+1 : 16);
		final HashSet<String> ids=new HashSet<String>(capacity);
		final ByteFile bf=ByteFile.makeByteFile(reps, true);
		boolean header=false, sequence=false;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length==0){continue;}
				if(line[0]=='>'){
					require(!header || sequence, "Empty representative FASTA record");
					int end=1;
					while(end<line.length && line[end]>32){end++;}
					require(end>1, "Empty representative FASTA identifier");
					String id=new String(line, 1, end-1, StandardCharsets.US_ASCII);
					if(!ids.add(id)){throw new IllegalArgumentException("Duplicate representative FASTA ID: "+id);}
					header=true;
					sequence=false;
				}else{
					require(header, "Sequence before first representative FASTA header");
					sequence=true;
				}
			}
		}finally{require(!bf.close(), "Representative reader reported an I/O failure");}
		require(header && sequence, "Representative FASTA is empty or ends with an empty record");
		log.println("Representative inventory: "+ids.size());
		return ids;
	}

	private void finishGroup(int size, int self){
		if(size<=0 || self!=1){throw new IllegalArgumentException("Each cluster must contain its representative once: size="+size+", self occurrences="+self);}
		if(size>=hist.length){hist=Arrays.copyOf(hist, Math.max(size+1, hist.length*2));}
		hist[size]++;
		families++;
		maxSize=Math.max(maxSize, size);
	}

	/** Writes descending size frequencies, with cumulative totals for whole tied bins. */
	private void writeReports(){
		long nf=0, ng=0;
		ByteStreamWriter bw=writer(".histogram.tsv");
		bw.print("cluster_size\tfamilies\tproteins\tcumulative_families\tcumulative_proteins\tgene_fraction\n");
		for(int size=maxSize; size>=1; size--){
			if(hist[size]==0){continue;}
			nf+=hist[size];
			ng+=hist[size]*size;
			bw.print(size).tab().print(hist[size]).tab().print(hist[size]*size).tab()
				.print(nf).tab().print(ng).tab().print((double)ng/rows, 9).nl();
		}
		close(bw);
		require(nf==families && ng==rows, "Histogram must conserve family and membership counts");
		bw=writer(".targets.tsv");
		bw.print("target_percent\tminimum_families\tcovered_proteins\tgene_fraction\tcutoff_cluster_size\n");
		for(int pct : TARGETS){
			final long target=(Math.multiplyExact(rows, pct)+99)/100;
			long count=0, genes=0;
			for(int size=maxSize; size>=1; size--){
				if(genes+hist[size]*size>=target){
					final long take=(target-genes+size-1)/size;
					count+=take;
					genes+=take*size;
					bw.print(pct).tab().print(count).tab().print(genes).tab().print((double)genes/rows, 9).tab().print(size).nl();
					break;
				}
				count+=hist[size];
				genes+=hist[size]*size;
			}
		}
		close(bw);
		bw=writer(".top.tsv");
		bw.print("families\tcovered_proteins\tgene_fraction\n");
		for(long wanted : TOP){
			if(wanted>families){continue;}
			long remaining=wanted, genes=0;
			for(int size=maxSize; size>0 && remaining>0; size--){
				final long take=Math.min(remaining, hist[size]);
				genes+=take*size;
				remaining-=take;
			}
			require(remaining==0, "Top-family request exceeded validated inventory");
			bw.print(wanted).tab().print(genes).tab().print((double)genes/rows, 9).nl();
		}
		close(bw);
		bw=writer(".bins.tsv");
		bw.print("minimum_size\tmaximum_size\tfamilies\tproteins\n");
		for(long lo=1, hi=1; lo<=maxSize; lo=hi+1, hi*=2){
			long fs=0, gs=0;
			for(int size=(int)lo; size<=Math.min(hi, maxSize); size++){fs+=hist[size]; gs+=hist[size]*size;}
			bw.print(lo).tab().print(hi).tab().print(fs).tab().print(gs).nl();
		}
		close(bw);
		bw=writer(".summary.tsv");
		bw.print("metric\tvalue\n").print("families\t").print(families).nl().print("proteins\t").print(rows).nl()
			.print("singletons\t").print(hist[1]).nl().print("largest_cluster\t").print(maxSize).nl();
		close(bw);
	}

	private ByteStreamWriter writer(String suffix){
		ByteStreamWriter bw=new ByteStreamWriter(out+suffix, false, false, true);
		bw.start();
		return bw;
	}

	private static void close(ByteStreamWriter bw){require(!bw.poisonAndWait(), "Output writer reported an I/O failure");}

	private static boolean same(byte[] row, int offset, int length, byte[] id){
		assert(offset>=0 && length>=0 && offset+length<=row.length) : "Identifier comparison must stay inside its parsed TSV field";
		if(id==null || length!=id.length){return false;}
		for(int i=0; i<length; i++){if(row[offset+i]!=id[i]){return false;}}
		return true;
	}

	private static void require(boolean ok, String message){if(!ok){throw new IllegalArgumentException(message);}}

	private String in, reps, out;
	private long expectedFamilies=-1, expectedGenes=-1, rows=0, families=0;
	private int maxSize=0;
	private long[] hist=new long[256];
	private final PrintStream log;
	private static final String[] SUFFIXES={".histogram.tsv", ".targets.tsv", ".top.tsv", ".bins.tsv", ".summary.tsv"};
	private static final int[] TARGETS={25, 50, 75, 80, 90, 95, 99, 100};
	private static final long[] TOP={1, 2, 3, 10, 100, 1000, 4000, 8000, 16000, 32000, 64000, 100000, 250000, 500000, 1000000};
}
