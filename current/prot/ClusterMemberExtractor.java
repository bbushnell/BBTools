package prot;

import java.io.File;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/**
 * Extracts selected clusters with memory bounded by one family's IDs, not the
 * entire selected corpus. MMseqs createtsv membership groups and result2flat
 * sequence groups may have different orders. The latter's empty header record
 * is a required cluster boundary, followed by a real representative record.
 * ByteFile is intentional: ordinary FASTA readers discard these empty markers.
 * Every emitted nonempty record consumes exactly one expected member ID; missing,
 * repeated or extra IDs fail. No sequence or membership cap is applied.
 * Raw intermediates and the seven-column manifest are written in a fresh root;
 * PASS appears only after all selected groups and optional source totals match.
 * @author Brian Bushnell, Keqing
 */
public final class ClusterMemberExtractor {

	public static void main(String[] args){
		ClusterMemberExtractor x=new ClusterMemberExtractor(args);
		try{x.process();}catch(Throwable t){t.printStackTrace(); System.exit(1);}
		finally{Shared.closeStream(x.log);}
	}

	private ClusterMemberExtractor(String[] args){
		PreParser pp=new PreParser(args, getClass(), false); log=pp.outstream;
		Parser parser=new Parser();
		for(String arg : pp.args){
			final int eq=arg.indexOf('=');
			final String a=(eq<0 ? arg : arg.substring(0, eq)).toLowerCase(Locale.ROOT);
			final String b=(eq<0 ? null : arg.substring(eq+1));
			if(a.equals("selected")){selected=b;}
			else if(a.equals("seqs")){seqs=b;}
			else if(a.equals("expectedgenes")){expectedGenes=Long.parseLong(b);}
			else if(a.equals("expectedfamilies")){expectedFamilies=Long.parseLong(b);}
			else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		in=parser.in1; out=parser.out1;
		require(in!=null && selected!=null && seqs!=null && out!=null && !parser.append,
			"Required: in=clusters.tsv selected=selection.tsv seqs=MMseqs_segmented.fasta out=new_directory; no append");
		require(!new File(out).exists(), "Output root already exists: "+out);
		Shared.setThreads(1); ByteFile.FORCE_MODE_BF1=true;
	}

	private void process() throws Exception{
		loadSelection();
		require(new File(out).mkdirs(), "Could not create fresh output directory: "+out);
		writeMemberships();
		extractSequences();
		final ByteStreamWriter manifest=writer("family_manifest.tsv");
		manifest.print("active_index\tfamily_id\trep_id\tassigned_count\tsequence_bytes\tfasta_file\tfasta_sha80\n");
		long total=0, residues=0;
		for(Family f : order){
			require(f.sequenceDone && f.actual==f.expected, "Selected sequence group missing/incomplete: "+f.rep);
			manifest.print(f.active).tab().print(f.id).tab().print(f.rep).tab().print(f.actual).tab().print(f.residues)
				.tab().print(f.fasta()).tab().print(DigestSuffix.file(out+"/"+f.fasta())).nl();
			total+=f.actual; residues+=f.residues;
		}
		close(manifest);
		ByteStreamWriter bw=writer("PASS");
		bw.print("CLUSTER_MEMBER_EXTRACTION_PASS families=").print(order.size()).print(" members=").print(total)
			.print(" residues=").print(residues).print(" source_records=").print(records).print(" source_groups=").print(boundaries).nl();
		close(bw);
		log.println("CLUSTER_MEMBER_EXTRACTION_PASS families="+order.size()+" members="+total+" residues="+residues);
	}

	private void loadSelection(){
		final ByteFile bf=ByteFile.makeByteFile(selected, true);
		final LineParser1 lp=new LineParser1('\t');
		final HashSet<String> ids=new HashSet<String>();
		boolean header=false;
		int previous=-1;
		try{
			for(byte[] row=bf.nextLine(); row!=null; row=bf.nextLine()){
				lp.set(row);
				if(!header){require(Arrays.equals(row, "active_index\tfamily_id\trep_id\tmember_count".getBytes(StandardCharsets.US_ASCII)), "Unsupported selection header"); header=true; continue;}
				require(lp.terms()==4, "Selection must have four columns");
				final Family f=new Family(lp.parseInt(0), lp.parseInt(1), lp.parseString(2), lp.parseInt(3));
				require(lp.termEquals(Integer.toString(f.id), 1), "Selected permanent IDs must use canonical decimal spelling: "+lp.parseString(1));
				require(f.active>=0 && (previous<0 || f.active==previous+1) && f.id>=0 && f.expected>0 && !f.rep.isEmpty(), "Invalid selection row: "+f.rep);
				require(ids.add(lp.parseString(1)) && byRep.put(f.rep, f)==null, "Duplicate selected family ID/representative: "+f.rep);
				order.add(f); previous=f.active;
			}
		}finally{require(!bf.close(), "Selection input I/O failure");}
		require(!order.isEmpty(), "Empty selection");
	}

	/** One open writer at a time; only selected cluster IDs are materialized on disk. */
	private void writeMemberships(){
		final ByteFile bf=ByteFile.makeByteFile(in, true);
		ByteStreamWriter bw=null;
		Family current=null;
		byte[] previous=null;
		long rows=0;
		int n=0;
		final ByteBuilder member=new ByteBuilder();
		try{
			for(byte[] row=bf.nextLine(); row!=null; row=bf.nextLine()){
				final int tab=AdditionalFamilySelector.separator(row, ++rows);
				if(!AdditionalFamilySelector.same(row, 0, tab, previous)){
					if(bw!=null){close(bw); bw=null; require(n==current.expected, "Membership count differs for "+current.rep+": "+n+" expected="+current.expected);}
					previous=Arrays.copyOf(row, tab);
					current=byRep.get(new String(previous, StandardCharsets.US_ASCII));
					n=0;
					if(current!=null){
						require(!current.membersDone, "Repeated selected membership group: "+current.rep);
						current.membersDone=true; bw=writer(current.ids());
					}
				}
				if(bw!=null){member.clear().append(row, tab+1, row.length-tab-1).nl(); bw.print(member); n++;}
			}
			if(bw!=null){close(bw); bw=null; require(n==current.expected, "Final membership count differs for "+current.rep);}
		}finally{
			if(bw!=null){close(bw);}
			require(!bf.close(), "Membership input I/O failure");
		}
		for(Family f : order){require(f.membersDone, "Selected cluster missing from membership input: "+f.rep);}
		require(expectedGenes<0 || rows==expectedGenes, "Membership source total differs: "+rows+" expected="+expectedGenes);
		log.println("Selected membership files complete; source rows="+rows);
	}

	/** Streams the special MMseqs format without retaining any nonselected sequences. */
	private void extractSequences(){
		final ByteFile bf=ByteFile.makeByteFile(seqs, true);
		byte[] header=null;
		boolean hasSequence=false, hasBoundary=false;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length==0){continue;}
				if(line[0]=='>'){
					validateHeader(line);
					if(header!=null && !hasSequence){
						finishGroup(); boundaries++;
						require(Arrays.equals(header, line), "Empty MMseqs boundary must be followed by its identical representative header");
						startGroup(new String(header, 1, header.length-1, StandardCharsets.US_ASCII));
						hasBoundary=true;
					}
					header=line; hasSequence=false;
				}else{
					require(header!=null && hasBoundary, "Sequence appeared before a valid MMseqs cluster boundary");
					if(!hasSequence){
						records++;
						if(current!=null){
							final String id=new String(header, 1, header.length-1, StandardCharsets.US_ASCII);
							require(remaining.remove(id), "Unexpected or repeated member in "+current.rep+": "+id);
							current.actual++; output.println(header);
						}
						if(records%10000000==0){log.println("Source proteins scanned: "+records+"; groups: "+boundaries);}
					}
					if(current!=null){
						for(byte b : line){if(!((b>='A' && b<='Z') || b=='*')){throw new IllegalArgumentException("Non-protein byte in selected family: "+current.rep);}}
						output.println(line); current.residues+=line.length;
					}
					hasSequence=true;
				}
			}
			require(header!=null && hasSequence, "MMseqs sequence input is empty or ends in an empty record");
			finishGroup();
		}finally{
			if(output!=null){close(output); output=null;}
			require(!bf.close(), "MMseqs sequence input I/O failure");
		}
		require(expectedGenes<0 || records==expectedGenes, "Sequence total differs: "+records+" expected="+expectedGenes);
		require(expectedFamilies<0 || boundaries==expectedFamilies, "Sequence group total differs: "+boundaries+" expected="+expectedFamilies);
	}

	private void startGroup(String rep){
		current=byRep.get(rep);
		if(current==null){return;}
		require(!current.sequenceDone, "Repeated selected sequence group: "+rep);
		remaining=new HashSet<String>();
		final ByteFile bf=ByteFile.makeByteFile(out+"/"+current.ids(), true);
		try{
			for(byte[] row=bf.nextLine(); row!=null; row=bf.nextLine()){
				String id=new String(row, StandardCharsets.US_ASCII);
				require(remaining.add(id), "Duplicate selected membership ID in "+rep+": "+id);
			}
		}finally{require(!bf.close(), "Selected membership input I/O failure");}
		require(remaining.size()==current.expected && remaining.contains(rep), "Selected membership must contain the representative and expected count: "+rep);
		output=writer(current.fasta());
	}

	private void finishGroup(){
		if(current==null){return;}
		assert(output!=null && remaining!=null) : "Selected group owns exactly one output writer and remaining-ID set";
		close(output); output=null;
		require(remaining.isEmpty() && current.actual==current.expected, "Missing member sequences in "+current.rep+": remaining="+remaining.size());
		current.sequenceDone=true; remaining=null; current=null;
	}

	private static void validateHeader(byte[] header){
		require(header.length>1, "Empty protein identifier");
		for(int i=1; i<header.length; i++){
			if(header[i]<=32 || header[i]>=127){throw new IllegalArgumentException("MMseqs headers must contain only one printable identifier");}
		}
	}
	private ByteStreamWriter writer(String name){ByteStreamWriter bw=new ByteStreamWriter(out+"/"+name, false, false, true); bw.start(); return bw;}
	private static void require(boolean ok, String message){AdditionalFamilySelector.require(ok, message);}
	private static void close(ByteStreamWriter bw){AdditionalFamilySelector.close(bw);}

	private static final class Family {
		Family(int a, int i, String r, int n){active=a; id=i; rep=r; expected=n;}
		String ids(){return "rank_"+active+".ids";}
		String fasta(){return "rank_"+active+".faa";}
		final int active, id, expected;
		final String rep;
		int actual;
		long residues;
		boolean membersDone, sequenceDone;
	}
	private String in, out, selected, seqs;
	private long expectedGenes=-1, expectedFamilies=-1, records=0, boundaries=0;
	private final PrintStream log;
	private final ArrayList<Family> order=new ArrayList<Family>();
	private final HashMap<String,Family> byRep=new HashMap<String,Family>();
	private Family current;
	private HashSet<String> remaining;
	private ByteStreamWriter output;
}
