package prot;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.atomic.AtomicInteger;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.Parse;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Builds one {@link AAGraph} consensus per pre-computed protein family (the hybrid
 * foundation step): keep an existing clustering's membership (e.g. mmseqs), but
 * replace each family's representative with the consensus of its members, so the
 * search representative fits divergent members instead of being one arbitrary
 * longest member. The consensus reps are then re-searched (by mmseqs) to rebuild
 * the feature cache.
 *
 * <p>Inputs: a {@code rep<TAB>member} cluster TSV (mmseqs
 * {@code *_cluster.tsv} format) and the member-sequence FASTA. Members are
 * subsampled per family ({@link #maxMembers}) and consensus-building runs in
 * parallel across families (each family's graph is independent). Output is a
 * FASTA whose headers are the family representative ids, so the downstream
 * family-to-column mapping is unchanged.</p>
 *
 * <p>Usage: {@code consensusrepbuilder.sh clusters=cl_cluster.tsv seqs=bac.faa
 * out=consensus_reps.fasta [families=familylist.tsv] [maxmembers=300] [t=64]}</p>
 *
 * @author Eru
 */
public final class ConsensusRepBuilder {

	/*--------------------------------------------------------------*/
	/*----------------            Config            ----------------*/
	/*--------------------------------------------------------------*/

	private String clustersFile, seqsFile, out, familiesFile;
	private boolean overwrite=true;
	/** Max members used to build a family's consensus (keep-first subsample). */
	private int maxMembers=300;
	/** Skip families with fewer than this many members. */
	private int minMembers=1;
	/** Skip members longer than this from a consensus (guards O(m*n) matrix RAM/OOM). */
	private int maxMemberLen=6000;
	private int threads=Shared.threads();
	/** Families whose consensus build threw (logged individually; distinct from skips). */
	private int crashedFamilies=0;

	/** AAGraph knobs (see {@link AAGraph}). */
	private int pad=20;
	private int consensusPasses=1;
	private boolean weightByIdentity=false;
	private float identityCeiling=40f;
	private float MAF_sub=0.25f, MAF_del=0.5f, MAF_ins=0.5f;
	/** Trim consensus ends below this fraction of max depth (strips smeared low-depth padding tails). */
	private float trimDepthFraction=0.1f;

	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		final Timer t=new Timer();
		final ConsensusRepBuilder x=new ConsensusRepBuilder(args);
		x.process(t);
	}

	public ConsensusRepBuilder(String[] args){
		if(args.length==0){printUsageAndExit();}
		for(String arg : args){
			final int eq=arg.indexOf('=');
			final String a=(eq<0 ? arg : arg.substring(0, eq)).toLowerCase();
			final String b=(eq<0 ? null : arg.substring(eq+1));
			if(a.equals("clusters") || a.equals("cluster") || a.equals("in")){clustersFile=b;}
			else if(a.equals("seqs") || a.equals("sequences") || a.equals("faa")){seqsFile=b;}
			else if(a.equals("out") || a.equals("o")){out=b;}
			else if(a.equals("families") || a.equals("familylist")){familiesFile=b;}
			else if(a.equals("maxmembers") || a.equals("max")){maxMembers=Integer.parseInt(b);}
			else if(a.equals("minmembers") || a.equals("min")){minMembers=Integer.parseInt(b);}
			else if(a.equals("maxmemberlen") || a.equals("maxlen")){maxMemberLen=Integer.parseInt(b);}
			else if(a.equals("t") || a.equals("threads")){threads=Integer.parseInt(b);}
			else if(a.equals("pad")){pad=Integer.parseInt(b);}
			else if(a.equals("passes")){consensusPasses=Integer.parseInt(b);}
			else if(a.equals("weightbyidentity") || a.equals("ani")){weightByIdentity=Parse.parseBoolean(b);}
			else if(a.equals("identityceiling") || a.equals("ceiling")){identityCeiling=Float.parseFloat(b);}
			else if(a.equals("maf_sub") || a.equals("mafsub")){MAF_sub=Float.parseFloat(b);}
			else if(a.equals("maf_del") || a.equals("mafdel")){MAF_del=Float.parseFloat(b);}
			else if(a.equals("maf_ins") || a.equals("mafins")){MAF_ins=Float.parseFloat(b);}
			else if(a.equals("trimdepth") || a.equals("trimdepthfraction")){trimDepthFraction=Float.parseFloat(b);}
			else if(a.equals("overwrite") || a.equals("ow")){overwrite=Parse.parseBoolean(b);}
			else if(a.equals("-h") || a.equals("--help") || a.equals("help")){printUsageAndExit();}
			else{throw new RuntimeException("Unknown argument: "+arg);}
		}
		if(clustersFile==null || seqsFile==null || out==null){
			System.err.println("Error: clusters=, seqs=, and out= are all required.\n");
			printUsageAndExit();
		}
		if(threads<1){threads=1;}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Pipeline           ----------------*/
	/*--------------------------------------------------------------*/

	private void process(final Timer t){
		final HashSet<String> keptReps=(familiesFile==null ? null : readKeptReps(familiesFile));
		if(keptReps!=null){System.err.println("Kept families: "+keptReps.size());}

		//Pass 1: group members per (kept) rep, subsampled; collect the needed member ids.
		final ArrayList<String> repOrder=new ArrayList<String>();
		final HashMap<String, ArrayList<String>> repMembers=new HashMap<String, ArrayList<String>>();
		final HashSet<String> needed=new HashSet<String>();
		readClusters(keptReps, repOrder, repMembers, needed);
		System.err.println("Families: "+repOrder.size()+"; distinct member seqs needed: "+needed.size());

		//Load only the needed member sequences (encoded), streaming the FASTA once.
		final HashMap<String, byte[]> encById=loadSequences(seqsFile, needed);
		System.err.println("Loaded "+encById.size()+" of "+needed.size()+" needed sequences.");

		//Build a consensus per family in parallel; results indexed by family order.
		final String[] consensus=buildAll(repOrder, repMembers, encById);

		writeFasta(repOrder, consensus, out, overwrite);

		t.stop();
		int written=0;
		for(String c : consensus){if(c!=null){written++;}}
		final int missing=repOrder.size()-written;
		System.err.println("Wrote "+written+" consensus reps ("+missing+" not written: "+
			crashedFamilies+" crashed, "+(missing-crashedFamilies)+" below min/no-members) in "+t);
		if(crashedFamilies>0){
			System.err.println("WARNING: "+crashedFamilies+" families crashed during consensus build "+
				"(see ERROR lines above); their reps are ABSENT from the output.");
		}
	}

	/** Header line prefix in a familylist TSV; checked without allocating the whole line. */
	private static final String FAMILYLIST_HEADER_PREFIX="rank\t";

	/** Reads the kept representative ids from a familylist TSV (column: cluster_rep). */
	private static HashSet<String> readKeptReps(final String file){
		final HashSet<String> set=new HashSet<String>();
		final ByteFile bf=ByteFile.makeByteFile(file, true);
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0 || line[0]=='#'){continue;}
			if(Tools.startsWith(line, FAMILYLIST_HEADER_PREFIX)){continue;}//header
			final int tab1=indexOf(line, (byte)'\t', 0);
			if(tab1<0){continue;}//fewer than 2 columns -- matches original parts.length>=2 guard
			final int tab2=indexOf(line, (byte)'\t', tab1+1);
			final int end=(tab2<0 ? line.length : tab2);
			set.add(new String(line, tab1+1, end-tab1-1));//column 1 (0-based): cluster_rep
		}
		bf.close();
		return set;
	}

	/**
	 * Streams the cluster TSV (rep&lt;TAB&gt;member), grouping subsampled members per kept rep.
	 * A real mmseqs cluster.tsv groups all of a representative's member rows consecutively, so
	 * this caches the previous line's rep bytes and resolved list/kept-status: a line whose rep
	 * field byte-matches the cached one reuses that list directly, skipping both the per-line
	 * {@code new String(rep)} allocation and the repMembers lookup that the original did on
	 * every single row (a family with 300 members otherwise re-allocates and re-looks-up the
	 * identical rep string 300 times). Falls back to the normal path whenever the cache misses
	 * (including the same rep reappearing non-consecutively), so correctness never depends on
	 * the input actually being grouped -- only the speed win does.
	 */
	private void readClusters(final HashSet<String> keptReps, final ArrayList<String> repOrder,
			final HashMap<String, ArrayList<String>> repMembers, final HashSet<String> needed){
		final ByteFile bf=ByteFile.makeByteFile(clustersFile, true);
		byte[] lastRepBytes=null;
		int lastRepLen=-1;
		ArrayList<String> lastList=null;
		boolean lastKept=false;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0 || line[0]=='#'){continue;}
			final int tab=indexOf(line, (byte)'\t', 0);
			if(tab<0){continue;}

			final ArrayList<String> list;
			if(lastRepBytes!=null && sameRegion(line, tab, lastRepBytes, lastRepLen)){
				if(!lastKept){continue;}
				list=lastList;
			}else{
				final String rep=new String(line, 0, tab);
				final boolean kept=(keptReps==null || keptReps.contains(rep));
				//Copy just the rep prefix (not a reference to the whole line) so the cache's
				//correctness never depends on which ByteFile implementation's nextLine() is
				//selected at runtime -- this runs once per distinct rep, not once per row.
				lastRepBytes=java.util.Arrays.copyOfRange(line, 0, tab); lastRepLen=tab; lastKept=kept;
				if(!kept){lastList=null; continue;}
				ArrayList<String> l=repMembers.get(rep);
				if(l==null){
					l=new ArrayList<String>();
					repMembers.put(rep, l);
					repOrder.add(rep);
				}
				list=l; lastList=l;
			}

			if(list.size()<maxMembers){//keep-first subsample (deterministic)
				final String member=new String(line, tab+1, line.length-tab-1);
				list.add(member);
				needed.add(member);
			}
		}
		bf.close();
	}

	/** True if {@code a[0,aLen)} equals {@code b[0,bLen)} byte-for-byte. */
	private static boolean sameRegion(final byte[] a, final int aLen, final byte[] b, final int bLen){
		if(aLen!=bLen){return false;}
		for(int i=0; i<aLen; i++){if(a[i]!=b[i]){return false;}}
		return true;
	}

	/** Streams the member FASTA, keeping only sequences whose id is needed (encoded leniently). */
	private static HashMap<String, byte[]> loadSequences(final String file, final HashSet<String> needed){
		final HashMap<String, byte[]> map=new HashMap<String, byte[]>(needed.size()*2);
		final FileFormat ff=FileFormat.testInput(file, FileFormat.FASTA, null, true, true);
		final ConcurrentReadInputStream cris=
				ConcurrentReadInputStream.getReadInputStream(-1, false, ff, null);
		cris.start();
		for(ListNum<Read> ln=cris.nextList(); ln!=null && ln.size()>0; ln=cris.nextList()){
			for(Read r : ln){
				if(r.bases==null){continue;}
				if(r.bases.length==0){continue;}//skip empty-sequence records (e.g. mmseqs separators)
				final String id=firstToken(r.id);
				if(needed.contains(id) && !map.containsKey(id)){
					map.put(id, encodeLenient(r.bases));
				}
			}
			cris.returnList(ln);
		}
		cris.close();
		return map;
	}

	/** Builds a consensus for every family, in parallel across families. */
	private String[] buildAll(final ArrayList<String> repOrder,
			final HashMap<String, ArrayList<String>> repMembers,
			final HashMap<String, byte[]> encById){
		final int n=repOrder.size();
		final String[] out=new String[n];
		final AtomicInteger cursor=new AtomicInteger(0);
		final AtomicInteger crashed=new AtomicInteger(0);
		final int nt=Math.min(threads, Math.max(1, n));
		final Thread[] pool=new Thread[nt];
		for(int i=0; i<nt; i++){
			pool[i]=new Thread(){
				@Override
				public void run(){
					for(int idx=cursor.getAndIncrement(); idx<n; idx=cursor.getAndIncrement()){
						try{
							out[idx]=buildOne(repMembers.get(repOrder.get(idx)), encById);
						}catch(Throwable e){
							//Log the specific family so a crash is never mistaken for a legitimate
							//low-membership skip (both would otherwise leave out[idx]==null).
							crashed.incrementAndGet();
							System.err.println("ERROR building family "+repOrder.get(idx)+": "+e);
						}
					}
				}
			};
			pool[i].start();
		}
		for(Thread th : pool){
			try{th.join();}catch(InterruptedException e){throw new RuntimeException(e);}
		}
		crashedFamilies=crashed.get();
		return out;
	}

	/** Builds one family's consensus string, or null if too few members resolved. */
	private String buildOne(final ArrayList<String> memberIds, final HashMap<String, byte[]> encById){
		final ArrayList<byte[]> members=new ArrayList<byte[]>();
		for(String id : memberIds){
			final byte[] e=encById.get(id);
			if(e!=null && e.length>0 && e.length<=maxMemberLen){members.add(e);}
		}
		if(members.size()<minMembers || members.isEmpty()){return null;}

		byte[] pivot=members.get(0);
		for(byte[] m : members){if(m.length>pivot.length){pivot=m;}}
		byte[] cons=graphConsensus(pivot, members);
		for(int p=1; p<consensusPasses && cons.length>0; p++){cons=graphConsensus(cons, members);}
		if(cons.length==0){cons=pivot;}
		return AAGraph.decode(cons);
	}

	/** Piles members onto an X-padded pivot graph and returns the trace. */
	private byte[] graphConsensus(final byte[] pivotEnc, final ArrayList<byte[]> members){
		final AAGraph g=new AAGraph(pivotEnc, pad);
		g.weightByIdentity=weightByIdentity;
		g.identityCeiling=identityCeiling;
		g.MAF_sub=MAF_sub; g.MAF_del=MAF_del; g.MAF_ins=MAF_ins;
		g.trimDepthFraction=trimDepthFraction;
		for(byte[] m : members){
			final AAAlignment aln=AAAligner.alignGlocal(m, g.pivot, true);
			if(aln!=null){g.add(m, aln);}
		}
		return g.traverse();
	}

	/** Writes the consensus reps as a FASTA keyed by family representative id. */
	private static void writeFasta(final ArrayList<String> repOrder, final String[] consensus,
			final String out, final boolean overwrite){
		final FileFormat ff=FileFormat.testOutput(out, FileFormat.FASTA, null, false, overwrite, false, false);
		final ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		final ByteBuilder bb=new ByteBuilder();
		for(int i=0; i<repOrder.size(); i++){
			if(consensus[i]==null){continue;}
			bb.clear();
			bb.append('>').append(repOrder.get(i)).nl();
			bb.append(consensus[i]).nl();
			bsw.print(bb.toBytes());
		}
		bsw.poisonAndWait();
	}

	/*--------------------------------------------------------------*/
	/*----------------           Helpers            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Encodes amino acids leniently for real gene-called proteins: standard
	 * residues map to 0-19; stops ('*') and gaps ('-','.') are dropped; any other
	 * letter (X/B/Z/J/U/O/...) maps to X. Whitespace and digits are ignored.
	 * @param raw Raw ASCII residues.
	 * @return Encoded residues (0-19 or X_CODE).
	 */
	static byte[] encodeLenient(final byte[] raw){
		final byte[] enc=new byte[raw.length];
		int n=0;
		for(final byte rb : raw){
			final int c=rb&0xFF;
			if(c=='*' || c=='-' || c=='.' || c<=' '){continue;}
			final byte std=(c<128 ? AminoAcid.acidToNumber[c] : -1);
			if(std>=0 && std<=19){enc[n++]=std;}
			else if((c>='A' && c<='Z') || (c>='a' && c<='z')){enc[n++]=Blosum62.X_CODE;}
			//else: digit or punctuation, ignored
		}
		return (n==raw.length ? enc : java.util.Arrays.copyOf(enc, n));
	}

	/** First whitespace-delimited token of a FASTA header. */
	static String firstToken(final String header){
		if(header==null){return "";}
		int i=0;
		while(i<header.length() && !Character.isWhitespace(header.charAt(i))){i++;}
		return header.substring(0, i);
	}

	/** Index of the first occurrence of a byte at or after fromIndex, or -1. */
	private static int indexOf(final byte[] a, final byte b, final int fromIndex){
		for(int i=fromIndex; i<a.length; i++){if(a[i]==b){return i;}}
		return -1;
	}

	private static void printUsageAndExit(){
		System.err.println(
			"ConsensusRepBuilder (BBTools prot) -- AAGraph consensus rep per pre-clustered family\n"+
			"\n"+
			"Usage: consensusrepbuilder.sh clusters=<rep_member.tsv> seqs=<members.faa> out=<reps.fasta>\n"+
			"\n"+
			"Required:\n"+
			"  clusters=  mmseqs-style cluster TSV (rep<TAB>member per line).\n"+
			"  seqs=      FASTA of member sequences (ids match the TSV members).\n"+
			"  out=       Output consensus FASTA (headers = family rep ids).\n"+
			"Optional:\n"+
			"  families=  familylist TSV; only its cluster_rep column's families are built.\n"+
			"  maxmembers=Members per family used for consensus (default 300, keep-first).\n"+
			"  minmembers=Skip families below this member count (default 1).\n"+
			"  maxlen=    Skip members longer than this from a consensus (default 6000; OOM guard).\n"+
			"  t=         Threads (default all).\n"+
			"  ani=       Identity-inverted weighting (t/f, default f).\n"+
			"  ceiling=   Identity ceiling for ani weighting (percent, default 40).\n"+
			"  passes=    Consensus refinement passes (default 1).\n"+
			"  trimdepth= Trim ends below this fraction of max depth (default 0.1; strips padding tails).\n"+
			"  pad=       Pivot X-padding each end (default 20).\n"+
			"  ow=        Overwrite output (t/f, default t).\n");
		System.exit(0);
	}
}
