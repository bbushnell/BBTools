package prot;

import java.io.File;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Locale;
import java.util.PriorityQueue;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.Parser;
import parse.PreParser;
import shared.Shared;

/**
 * Selects the largest unrepresented original MMseqs clusters. Current permanent
 * family IDs and their original parent IDs resolve through the original ranked
 * family list; split children therefore exclude their original ancestor too.
 * Ranking is descending membership count, then ascending representative ID.
 * New active indexes follow the current roster; new permanent IDs follow its
 * maximum permanent ID. Output is experimental and does not mutate that roster.
 * The TSV must contain one contiguous group per representative, including its
 * self-member exactly once. Counts describe rows, not distinct biological genes.
 * @author Brian Bushnell, Keqing
 */
public final class AdditionalFamilySelector {

	public static void main(String[] args){
		AdditionalFamilySelector x=new AdditionalFamilySelector(args);
		try{x.process();}finally{Shared.closeStream(x.log);}
	}

	private AdditionalFamilySelector(String[] args){
		PreParser pp=new PreParser(args, getClass(), false);
		log=pp.outstream;
		Parser parser=new Parser();
		for(String arg : pp.args){
			final int eq=arg.indexOf('=');
			final String a=(eq<0 ? arg : arg.substring(0, eq)).toLowerCase(Locale.ROOT);
			final String b=(eq<0 ? null : arg.substring(eq+1));
			if(a.equals("original")){original=b;}
			else if(a.equals("roster")){roster=b;}
			else if(a.equals("count")){count=Integer.parseInt(b);}
			else if(a.equals("expectedfamilies")){expectedFamilies=Long.parseLong(b);}
			else if(a.equals("expectedgenes")){expectedGenes=Long.parseLong(b);}
			else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		in=parser.in1; out=parser.out1;
		require(in!=null && original!=null && roster!=null && out!=null && count>0,
			"Required: in=clusters.tsv original=original_familylist.tsv roster=current_roster.tsv out=prefix count=N>0");
		require(!parser.append, "Appending would mix independently ranked selections");
		for(String suffix : new String[]{".selected.tsv", ".excluded.tsv", ".summary.tsv"}){
			require(!new File(out+suffix).exists(), "Refusing existing output: "+out+suffix);
		}
		Shared.setThreads(1); ByteFile.FORCE_MODE_BF1=true;
	}

	private void process(){
		final ArrayList<String> roots=loadOriginal();
		final boolean[] represented=loadRoster(roots.size());
		final HashSet<String> excluded=new HashSet<String>();
		for(int i=0; i<roots.size(); i++){if(represented[i]){excluded.add(roots.get(i));}}
		final HashSet<String> seen=new HashSet<String>();
		final PriorityQueue<Family> largest=new PriorityQueue<Family>(count, BEST.reversed());
		final ByteFile bf=ByteFile.makeByteFile(in, true);
		byte[] previous=null;
		String rep=null;
		int n=0, self=0;
		try{
			for(byte[] row=bf.nextLine(); row!=null; row=bf.nextLine()){
				rows++;
				final int tab=separator(row, rows);
				if(!same(row, 0, tab, previous)){
					if(rep!=null){accept(rep, n, self, excluded, largest);}
					previous=Arrays.copyOf(row, tab);
					rep=new String(previous, StandardCharsets.US_ASCII);
					require(seen.add(rep), "Noncontiguous repeated cluster: "+rep);
					n=0; self=0;
				}
				n=Math.addExact(n, 1);
				if(same(row, tab+1, row.length-tab-1, previous)){self++;}
				if(rows%10000000==0){log.println("Membership rows: "+rows+"; groups: "+seen.size());}
			}
		}finally{require(!bf.close(), "Membership input I/O failure");}
		if(rep!=null){accept(rep, n, self, excluded, largest);}
		require(rows>0 && largest.size()==count, "Not enough unrepresented clusters: selected="+largest.size()+", requested="+count);
		require(expectedFamilies<0 || seen.size()==expectedFamilies, "Cluster total differs: "+seen.size()+" expected="+expectedFamilies);
		require(expectedGenes<0 || rows==expectedGenes, "Membership total differs: "+rows+" expected="+expectedGenes);
		for(String id : roots){require(seen.contains(id), "Original selected representative absent from cluster corpus: "+id);}
		ArrayList<Family> selected=new ArrayList<Family>(largest);
		Collections.sort(selected, BEST);
		ByteStreamWriter bw=writer(".selected.tsv");
		bw.print("active_index\tfamily_id\trep_id\tmember_count\n");
		long selectedMembers=0;
		for(int i=0; i<selected.size(); i++){
			Family f=selected.get(i); selectedMembers+=f.members;
			bw.print(Math.addExact(activeCount, i)).tab().print(Math.addExact(maxId+1, i)).tab().print(f.rep).tab().print(f.members).nl();
		}
		close(bw);
		bw=writer(".excluded.tsv"); bw.print("original_family_id\toriginal_rep_id\n");
		for(int i=0; i<roots.size(); i++){if(represented[i]){bw.print(i).tab().print(roots.get(i)).nl();}}
		close(bw);
		bw=writer(".summary.tsv");
		bw.print("metric\tvalue\ncurrent_families\t").print(activeCount).print("\nexcluded_original_roots\t").print(excluded.size())
			.print("\nsource_families\t").print(seen.size()).print("\nsource_memberships\t").print(rows)
			.print("\nselected_families\t").print(selected.size()).print("\nselected_memberships\t").print(selectedMembers)
			.print("\nsmallest_selected_family\t").print(selected.get(selected.size()-1).members).nl();
		close(bw);
		log.println("ADDITIONAL_FAMILY_SELECTION_PASS selected="+selected.size()+" members="+selectedMembers+" excluded_roots="+excluded.size());
	}

	/** The original list defines permanent IDs by contiguous rank, independently of the current roster. */
	private ArrayList<String> loadOriginal(){
		final ArrayList<String> roots=new ArrayList<String>();
		final HashSet<String> ids=new HashSet<String>();
		final LineParser1 lp=new LineParser1('\t');
		final ByteFile bf=ByteFile.makeByteFile(original, true);
		try{
			for(byte[] row=bf.nextLine(); row!=null; row=bf.nextLine()){
				if(row.length==0 || row[0]=='#'){continue;}
				lp.set(row);
				require(lp.terms()==3 && lp.parseInt(0)==roots.size(), "Original family list must have rank,rep_id,occ_total in contiguous rank order");
				final String id=lp.parseString(1);
				require(!id.isEmpty() && ids.add(id), "Empty or duplicate original representative: "+id);
				roots.add(id);
			}
		}finally{require(!bf.close(), "Original family list I/O failure");}
		require(!roots.isEmpty(), "Empty original family list");
		return roots;
	}

	/** Supports original IDs or direct original-parent IDs; unknown lineage fails instead of guessing. */
	private boolean[] loadRoster(int rootCount){
		final boolean[] represented=new boolean[rootCount];
		final HashSet<String> ids=new HashSet<String>(), reps=new HashSet<String>();
		final LineParser1 lp=new LineParser1('\t');
		final ByteFile bf=ByteFile.makeByteFile(roster, true);
		boolean header=false;
		try{
			for(byte[] row=bf.nextLine(); row!=null; row=bf.nextLine()){
				if(row.length==0 || row[0]=='#'){continue;}
				lp.set(row);
				if(!header){
					require(lp.terms()==18 && lp.termEquals("active_index", 0) && lp.termEquals("family_id", 1)
						&& lp.termEquals("rep_id", 2) && lp.termEquals("parent_family_id", 4), "Unsupported roster header");
					header=true; continue;
				}
				require(lp.terms()==18 && lp.parseInt(0)==activeCount, "Roster active indexes must be contiguous and ordered");
				final int id=lp.parseInt(1), root=lp.termEquals("NA", 4) ? id : lp.parseInt(4);
				require(lp.termEquals(Integer.toString(id), 1) && (lp.termEquals("NA", 4) || lp.termEquals(Integer.toString(root), 4)),
					"Roster IDs must use canonical decimal spelling at active index "+activeCount);
				require(id>=0 && id<Integer.MAX_VALUE-count && ids.add(lp.parseString(1)) && reps.add(lp.parseString(2)),
					"Invalid or duplicate permanent ID/representative at active index "+activeCount);
				require(root>=0 && root<rootCount, "Unresolved original ancestor: family="+id+", parent="+root);
				represented[root]=true; maxId=Math.max(maxId, id); activeCount++;
			}
		}finally{require(!bf.close(), "Current roster I/O failure");}
		require(header && activeCount>0, "Empty current roster");
		return represented;
	}

	private void accept(String rep, int n, int self, HashSet<String> excluded, PriorityQueue<Family> largest){
		require(n>0 && self==1, "Cluster must have exactly one representative member: "+rep+" size="+n+" self="+self);
		if(excluded.contains(rep)){return;}
		Family f=new Family(rep, n);
		if(largest.size()<count){largest.add(f);}
		else if(BEST.compare(f, largest.peek())<0){largest.poll(); largest.add(f);}
	}

	/** Strict two-field ID input shared with the selected-member extractor. */
	static int separator(byte[] row, long line){
		int tab=-1;
		for(int i=0; i<row.length; i++){
			if(row[i]=='\t'){
				if(tab>=0){throw new IllegalArgumentException("Extra membership field at row "+line);}
				tab=i;
			}else if(row[i]<=32 || row[i]>=127){throw new IllegalArgumentException("Invalid membership ID byte at row "+line);}
		}
		if(tab<=0 || tab>=row.length-1){throw new IllegalArgumentException("Expected rep<TAB>member at row "+line);}
		return tab;
	}

	static boolean same(byte[] row, int offset, int length, byte[] id){
		assert(offset>=0 && length>=0 && offset+length<=row.length) : "ID comparison must remain inside the parsed input field";
		if(id==null || length!=id.length){return false;}
		for(int i=0; i<length; i++){if(row[offset+i]!=id[i]){return false;}}
		return true;
	}

	private ByteStreamWriter writer(String suffix){ByteStreamWriter bw=new ByteStreamWriter(out+suffix, false, false, true); bw.start(); return bw;}
	static void close(ByteStreamWriter bw){require(!bw.poisonAndWait(), "Output writer I/O failure");}
	static void require(boolean ok, String message){if(!ok){throw new IllegalArgumentException(message);}}

	private static final class Family {
		Family(String r, int n){rep=r; members=n;}
		final String rep;
		final int members;
	}
	private static final Comparator<Family> BEST=new Comparator<Family>(){
		@Override public int compare(Family a, Family b){int n=Integer.compare(b.members, a.members); return n==0 ? a.rep.compareTo(b.rep) : n;}
	};
	private String original, roster, in, out;
	private int count=4000, activeCount=0, maxId=-1;
	private long rows, expectedFamilies=-1, expectedGenes=-1;
	private final PrintStream log;
}
