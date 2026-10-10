package prot;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Compares pinned old and reassigned family memberships by query ID and encoded
 * sequence, rather than count or FASTA ordering. Reads every declared record
 * strictly, including zero-member new families; malformed proteins never skip.
 * EMPTY is a reported state, not a decision to delete the old model.
 * @author Keqing
 */
public final class HbmMembershipCompare {

	public static void main(String[] args){
		try{
			Shared.setThreads(1); Shared.AMINO_IN=true; Read.VALIDATE_IN_CONSTRUCTOR=false;
			final PreParser pp=new PreParser(args, HbmMembershipCompare.class, false);
			try{run(HmmComparisonData.options(pp.args));}finally{Shared.closeStream(pp.outstream);}
		}catch(Throwable failure){failure.printStackTrace(); System.exit(1);}
	}

	private static void run(HashMap<String,String> o) throws Exception{
		for(String key : o.keySet()){
			require(Arrays.asList("old", "oldsha80", "new", "newsha80", "ranks", "out").contains(key), "Unknown membership comparison parameter");
		}
		final String oldFile=HmmComparisonData.required(o, "old"), newFile=HmmComparisonData.required(o, "new");
		HbmProfilePilot.requireHash(oldFile, HmmComparisonData.required(o, "oldsha80"));
		HbmProfilePilot.requireHash(newFile, HmmComparisonData.required(o, "newsha80"));
		final List<String[]> oldRows=HmmComparisonData.rows(oldFile), newRows=HmmComparisonData.rows(newFile);
		require(oldRows.size()>1 && oldRows.size()==newRows.size(), "Old/new family roster sizes differ");
		require(String.join("\t", oldRows.get(0)).equals(OLD_HEADER) && String.join("\t", newRows.get(0)).equals(NEW_HEADER), "Invalid member manifest schemas");
		final HashSet<String> ids=new HashSet<String>(), reps=new HashSet<String>();
		for(int i=1; i<oldRows.size(); i++){
			final String[] a=oldRows.get(i), b=newRows.get(i); final String rank=Integer.toString(i-1);
			require(a.length==7 && b.length==9 && a[0].equals(rank) && b[0].equals(rank) && a[1].equals(b[1]) && a[2].equals(b[2]), "Member manifests have different family identities");
			final int id=Integer.parseInt(a[1]);
			require(id>=0 && a[1].equals(Integer.toString(id)) && ids.add(a[1]) && reps.add(a[2]), "Invalid or repeated permanent family ID/representative");
			final long oldCount=Long.parseLong(a[3]), count=Long.parseLong(b[3]), raw=Long.parseLong(b[4]), encoded=Long.parseLong(b[5]);
			require(oldCount>0 && count>=0 && raw>=encoded && encoded>=count && (count>0 || raw==0) && b[8].equals(count==0 ? "EMPTY" : "NONEMPTY"), "Invalid old/new family membership totals/state");
			require(Paths.get(a[5]).isAbsolute() && Paths.get(b[6]).isAbsolute(), "Membership inputs need absolute pinned paths");
		}
		final int[] ranks=HbmProfileLibrary.readRanks(o.get("ranks"), oldRows.size()-1); Arrays.sort(ranks);
		final Path out=Paths.get(HmmComparisonData.required(o, "out"));
		require(!Files.exists(out), "Membership comparison output must be fresh"); Files.createDirectory(out);
		final ByteStreamWriter table=HmmComparisonData.writer(out.resolve("changes.tsv").toString());
		table.println("active_index\tfamily_id\trep_id\told_members\tnew_members\tadded_ids\tremoved_ids\tchanged_sequences\tstate");
		final ByteBuilder row=new ByteBuilder(); int changed=0, unchanged=0, empty=0; long added=0, removed=0, sequences=0;
		try{
			for(int rank : ranks){
				final String[] a=oldRows.get(rank+1), b=newRows.get(rank+1);
				final List<ProteinSequence> oldMembers=readMembers(a[5], a[6], Long.parseLong(a[3]), Long.parseLong(a[4]), -1);
				final List<ProteinSequence> newMembers=readMembers(b[6], b[7], Long.parseLong(b[3]), Long.parseLong(b[4]), Long.parseLong(b[5]));
				final long[] delta=compare(oldMembers, newMembers);
				final String state=newMembers.isEmpty() ? "EMPTY" : delta[0]+delta[1]+delta[2]==0 ? "UNCHANGED" : "CHANGED";
				if(state.equals("EMPTY")){empty++;}else if(state.equals("UNCHANGED")){unchanged++;}else{changed++;}
				added+=delta[0]; removed+=delta[1]; sequences+=delta[2];
				table.println(row.clear().append(rank).tab().append(a[1]).tab().append(a[2]).tab().append(oldMembers.size()).tab().append(newMembers.size())
					.tab().append(delta[0]).tab().append(delta[1]).tab().append(delta[2]).tab().append(state));
				System.err.println("HBM_MEMBERSHIP_COMPARE_PROGRESS rank="+rank+" state="+state);
			}
		}finally{HmmComparisonData.close(table);}
		HbmProfilePilot.requireHash(oldFile, o.get("oldsha80")); HbmProfilePilot.requireHash(newFile, o.get("newsha80"));
		final ByteStreamWriter summary=HmmComparisonData.writer(out.resolve("summary.tsv").toString());
		summary.println("families\tunchanged\tchanged_nonempty\tempty\tadded_ids\tremoved_ids\tchanged_sequences");
		summary.println(row.clear().append(ranks.length).tab().append(unchanged).tab().append(changed).tab().append(empty).tab().append(added).tab().append(removed).tab().append(sequences));
		HmmComparisonData.close(summary);
		final ByteStreamWriter recipe=HmmComparisonData.writer(out.resolve("recipe.tsv").toString());
		recipe.println("old_manifest_sha80\t"+o.get("oldsha80")+"\nnew_manifest_sha80\t"+o.get("newsha80")+"\nequality\texact_query_id_and_encoded_sequence");
		HmmComparisonData.close(recipe);
		final ByteStreamWriter pass=HmmComparisonData.writer(out.resolve("PASS").toString()); pass.println("HBM_MEMBERSHIP_COMPARE_PASS"); HmmComparisonData.close(pass);
	}

	/** Reads strictly and checks raw as well as encoded totals; an empty new family is valid. */
	static List<ProteinSequence> readMembers(String file, String pin, long expected, long rawExpected, long encodedExpected){
		require(expected>=0 && expected<=Integer.MAX_VALUE && rawExpected>=0, "Unsupported family dimensions");
		HbmProfilePilot.requireHash(file, pin);
		final ArrayList<ProteinSequence> members=new ArrayList<ProteinSequence>();
		final Streamer input=StreamerFactory.makeStreamer(FileFormat.testInput(file, FileFormat.FASTA, null, true, false), null, true, -1);
		long raw=0, encoded=0; input.start();
		try{
			for(ListNum<Read> list=input.nextList(); list!=null && list.size()>0; list=input.nextList()){
				for(Read read : list){
					require(read.mate==null && read.bases!=null && read.bases.length>0, "Empty or paired membership record");
					final ProteinSequence p=new ProteinSequence(HbmCompetitiveAssign.queryId(read.id), read.bases);
					members.add(p); raw+=read.bases.length; encoded+=p.length();
				}
			}
		}finally{input.close(); require(!input.errorState(), "Member FASTA reader failed");}
		require(members.size()==expected && raw==rawExpected && (encodedExpected<0 || encoded==encodedExpected), "Member FASTA count/raw/encoded totals differ");
		HbmProfilePilot.requireHash(file, pin); return members;
	}

	/** Sorts caller-owned lists and returns added IDs, removed IDs, and changed shared sequences. */
	static long[] compare(List<ProteinSequence> oldMembers, List<ProteinSequence> newMembers){
		checkUnique(oldMembers); checkUnique(newMembers);
		int a=0, b=0; final long[] out=new long[3];
		while(a<oldMembers.size() && b<newMembers.size()){
			final ProteinSequence x=oldMembers.get(a), y=newMembers.get(b); final int cmp=x.id.compareTo(y.id);
			if(cmp<0){out[1]++; a++;}else if(cmp>0){out[0]++; b++;}
			else{if(!Arrays.equals(x.enc, y.enc)){out[2]++;}a++; b++;}
		}
		out[1]+=oldMembers.size()-a; out[0]+=newMembers.size()-b;
		assert(oldMembers.size()+out[0]-out[1]==newMembers.size()) : "Membership ID merge must conserve old plus additions minus removals";
		return out;
	}
	private static void checkUnique(List<ProteinSequence> members){
		members.sort(Comparator.comparing(p->p.id));
		for(int i=1; i<members.size(); i++){
			if(members.get(i-1).id.equals(members.get(i).id)){throw new IllegalArgumentException("Duplicate member query ID: "+members.get(i).id);}
		}
	}
	private static void require(boolean ok, String reason){if(!ok){throw new IllegalArgumentException(reason);}}
	private static final String OLD_HEADER="active_index\tfamily_id\trep_id\tassigned_count\tsequence_bytes\tfasta_file\tfasta_sha80";
	private static final String NEW_HEADER="active_index\tfamily_id\trep_id\tassigned_count\tderived_raw_residues\tencoded_residues\tfasta_file\tfasta_sha80\tmembership_status";
}
