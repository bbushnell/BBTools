package prot;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import fileIO.ByteFile;
import parse.LineParser1;

/**
 * Independent, external verification of one IdentityGroupBuilder output pair (cluster.tsv +
 * sidecar.tsv) BEFORE the source edge file is deleted (Elly's review, 2026-08-31 -- three
 * rounds). Never trusts IdentityGroupBuilder's own exit code alone, and does NOT merely
 * re-derive aggregate statistics from the sidecar's own claims -- every check here is computed
 * from scratch against cluster.tsv's actual bytes and rows, then compared to the sidecar.
 * <p>
 * <b>Round 1 gaps found by review and fixed here:</b>
 * <ol>
 * <li>The whole-partition hash is now a RAW byte-stream digest of cluster.tsv (via a plain
 * FileInputStream, byte-for-byte, including any blank/malformed line) -- the first draft
 * reconstructed lines via ByteFile (which skips blank lines) before hashing, so a corruption
 * that inserted or removed a blank line would NOT have changed the computed hash. This is the
 * one check that must see literally every byte, so it deliberately does NOT use the same
 * line-based parsing as the rest of this class.</li>
 * <li>Per-group membership is now independently RECONSTRUCTED from cluster.tsv (grouping every
 * member row by its rep column) and each group's size + members_sha256 is recomputed using the
 * IDENTICAL algorithm IdentityGroupBuilder.writeSidecar uses (sorted member list, each member's
 * bytes + a newline byte fed to the digest in order) -- then every sidecar row must match a
 * reconstructed group EXACTLY on (rep, size, hash), with no missing, extra, or duplicate rep in
 * either direction. The first draft only checked the SUM of sizes and the group-ROW COUNT, so
 * two groups with swapped/wrong per-group size or hash entries (same totals, wrong assignment)
 * would have silently passed.</li>
 * </ol>
 * Also rejects a duplicate member_id in cluster.tsv (same member appearing in two rows, whether
 * under the same or a different rep) and any malformed (not-exactly-2-field) row.
 * <p>
 * <b>Round 3 refactor (Elly's review, 2026-08-31 -- corpus-wide gate):</b> the core verification
 * logic is now {@link #verify(String, String, int)}, a reusable library method returning a
 * {@link Result} rather than calling {@code System.exit} directly, so a corpus-wide caller (one
 * process verifying all 4432 families in one pass, e.g. a completeness/conservation gate) can
 * call the SAME logic this class's own {@code main()} uses per-family, instead of duplicating
 * it. {@code main()} is now a thin wrapper: call {@link #verify}, print, exit accordingly.
 * <p>
 * Exits 0 (silent on stdout beyond the PASS line) on PASS; prints one FAIL line per violated
 * check and exits 1 otherwise (every check always runs and is reported, never short-circuits on
 * the first failure) so a batch job's slurm log shows the complete picture.
 *
 * <p>Usage: {@code java -ea prot.IdentityGroupOutputVerifier cluster=<rank>.cluster.tsv
 *        sidecar=<rank>.cluster.tsv.sidecar.tsv expectedmembers=<N>}
 *
 * @author Eru
 */
public class IdentityGroupOutputVerifier {

	public static void main(String[] args){
		String clusterFile=null, sidecarFile=null;
		int expectedMembers=-1;
		for(String arg : args){
			final int eq=arg.indexOf('=');
			final String a=(eq<0 ? arg : arg.substring(0, eq)).toLowerCase();
			final String b=(eq<0 ? null : arg.substring(eq+1));
			if(a.equals("cluster")){clusterFile=b;}
			else if(a.equals("sidecar")){sidecarFile=b;}
			else if(a.equals("expectedmembers")){expectedMembers=Integer.parseInt(b);}
			else{throw new RuntimeException("Unknown argument: "+arg);}
		}
		if(clusterFile==null || sidecarFile==null || expectedMembers<0){
			throw new RuntimeException("Required: cluster=<cluster.tsv> sidecar=<sidecar.tsv> "
				+"expectedmembers=<N>");
		}

		final Result r=verify(clusterFile, sidecarFile, expectedMembers);
		if(!r.pass){
			for(String f : r.failures){System.err.println("VERIFY_FAIL: "+f);}
			System.err.println("IdentityGroupOutputVerifier: FAIL ("+r.failures.size()+" check(s) failed) for "
				+clusterFile+" / "+sidecarFile);
			System.exit(1);
		}
		System.err.println("IdentityGroupOutputVerifier: PASS -- "+r.clusterRows+" rows, "
			+r.groupCount+" groups, sha256="+r.rawHash);
	}

	/** Everything the per-family verification proves, for a caller that needs more than a bare
	 * pass/fail -- e.g. a corpus-wide gate that also wants group_count/largest_group_size to
	 * compute the fragmentation proxy in the same pass, without re-reading cluster.tsv again. */
	public static final class Result {
		public boolean pass;
		public final ArrayList<String> failures=new ArrayList<String>();
		public int clusterRows;
		public int groupCount;
		public int largestGroupSize;
		public String rawHash;
	}

	/** Runs the FULL independent verification of one cluster.tsv+sidecar pair and returns every
	 * result a caller might need -- never exits, never prints; the caller decides what to do
	 * with a {@link Result}. This is the SAME logic {@code main()} uses per-family; extracted so
	 * a corpus-wide caller (e.g. a completeness/conservation gate over all 4432 families) calls
	 * this once per family instead of duplicating the reconstruction/comparison logic. */
	public static Result verify(String clusterFile, String sidecarFile, int expectedMembers){
		final Result result=new Result();
		final ArrayList<String> failures=result.failures;

		//--- Check 1: RAW byte-stream hash of cluster.tsv (every byte, no line reconstruction). ---
		final String rawHash=hashFileRaw(clusterFile);
		result.rawHash=rawHash;

		//--- Check 2: parse cluster.tsv into rep->members groups, rejecting malformed/duplicate rows. ---
		final HashMap<String,ArrayList<String>> groups=new HashMap<String,ArrayList<String>>();
		final HashSet<String> seenMembers=new HashSet<String>();
		int clusterRows=0;
		{
			final ByteFile bf=ByteFile.makeByteFile(clusterFile, true);
			final LineParser1 lp=new LineParser1((byte)'\t');
			long lineNo=0;
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lineNo++;
				if(line.length==0){
					failures.add("cluster.tsv line "+lineNo+" is blank.");
					continue;
				}
				lp.set(line);
				if(lp.terms()!=2){
					failures.add("cluster.tsv line "+lineNo+" does not have exactly 2 fields (rep,member): "+new String(line));
					continue;
				}
				final String rep=lp.parseString(0), member=lp.parseString(1);
				if(!seenMembers.add(member)){
					failures.add("cluster.tsv has duplicate member '"+member+"' (line "+lineNo+") -- "
						+"membership rows must be unique.");
				}
				ArrayList<String> list=groups.get(rep);
				if(list==null){list=new ArrayList<String>(); groups.put(rep, list);}
				list.add(member);
				clusterRows++;
			}
			if(bf.close()){failures.add("ByteFile reported an I/O error reading "+clusterFile);}
		}
		result.clusterRows=clusterRows;
		if(clusterRows!=expectedMembers){
			failures.add("cluster.tsv has "+clusterRows+" valid rows, expected exactly "+expectedMembers
				+" (one per member).");
		}

		//--- Recompute each group's size + members_sha256 EXACTLY as IdentityGroupBuilder.writeSidecar does. ---
		final HashMap<String,Integer> recomputedSize=new HashMap<String,Integer>();
		final HashMap<String,String> recomputedHash=new HashMap<String,String>();
		int largestGroupSize=0;
		for(java.util.Map.Entry<String,ArrayList<String>> e : groups.entrySet()){
			final ArrayList<String> members=e.getValue();
			Collections.sort(members);//IdentityGroupBuilder sorts each group's members before hashing/writing
			final MessageDigest gd=sha256();
			for(String m : members){gd.update(m.getBytes()); gd.update((byte)'\n');}
			recomputedSize.put(e.getKey(), members.size());
			recomputedHash.put(e.getKey(), toHex(gd.digest()));
			if(members.size()>largestGroupSize){largestGroupSize=members.size();}
		}
		result.groupCount=groups.size();
		result.largestGroupSize=largestGroupSize;

		//--- Check 3: sidecar schema + exact per-group (rep,size,hash) match, both directions. ---
		String sidecarWholeHash=null;
		final HashSet<String> sidecarReps=new HashSet<String>();
		boolean sawHeader2=false;
		{
			final ByteFile bf=ByteFile.makeByteFile(sidecarFile, true);
			int lineNo=0;
			for(byte[] lineBytes=bf.nextLine(); lineBytes!=null; lineBytes=bf.nextLine()){
				lineNo++;
				final String line=new String(lineBytes);
				if(lineNo==1){
					final String[] f=line.split("\t");
					if(f.length!=2 || !f[0].equals("#whole_partition_sha256")){
						failures.add("sidecar line 1 malformed (expected '#whole_partition_sha256\\t<hash>'): "+line);
					}else{
						sidecarWholeHash=f[1];
					}
				}else if(lineNo==2){
					if(!line.equals("#rep\tsize\tmembers_sha256")){
						failures.add("sidecar line 2 does not match the expected header exactly: '"+line+"'");
					}else{
						sawHeader2=true;
					}
				}else{
					final String[] f=line.split("\t");
					if(f.length!=3){
						failures.add("sidecar group row "+lineNo+" does not have exactly 3 fields: "+line);
						continue;
					}
					final String rep=f[0];
					if(!sidecarReps.add(rep)){
						failures.add("sidecar has a DUPLICATE rep '"+rep+"' (row "+lineNo+").");
						continue;
					}
					final int size;
					try{size=Integer.parseInt(f[1]);}
					catch(NumberFormatException ex){
						failures.add("sidecar group row "+lineNo+" has a non-integer size field: "+line);
						continue;
					}
					if(size<=0){
						failures.add("sidecar group row "+lineNo+" has a nonpositive size ("+size+"): "+line);
					}
					if(!f[2].matches("[0-9a-f]{64}")){
						failures.add("sidecar group row "+lineNo+" members_sha256 is not exactly 64 lowercase hex chars: '"+f[2]+"'");
					}
					//The load-bearing check (Elly's round-2 review): this sidecar row must match an
					//INDEPENDENTLY RECOMPUTED group on (rep,size,hash) -- not just be internally
					//well-formed. A swapped/wrong per-group entry with correct aggregate totals is
					//exactly what a sum/count-only check (round 1) would have missed.
					if(!recomputedSize.containsKey(rep)){
						failures.add("sidecar rep '"+rep+"' (row "+lineNo+") does not exist as a "
							+"reconstructed group in cluster.tsv at all.");
					}else{
						final int trueSize=recomputedSize.get(rep);
						final String trueHash=recomputedHash.get(rep);
						if(size!=trueSize){
							failures.add("sidecar rep '"+rep+"' claims size "+size
								+" but the reconstructed group from cluster.tsv has size "+trueSize+".");
						}
						if(!f[2].equals(trueHash)){
							failures.add("sidecar rep '"+rep+"' claims members_sha256='"+f[2]
								+"' but the reconstructed group's true hash is '"+trueHash+"'.");
						}
					}
				}
			}
			if(bf.close()){failures.add("ByteFile reported an I/O error reading "+sidecarFile);}
			if(lineNo<2){failures.add("sidecar has fewer than 2 lines (missing header rows entirely).");}
		}
		if(!sawHeader2){/* already reported above */}
		if(sidecarWholeHash!=null && !sidecarWholeHash.equals(rawHash)){
			failures.add("sidecar whole_partition_sha256 ("+sidecarWholeHash+") does NOT match the "
				+"RAW byte-stream sha256 of cluster.tsv ("+rawHash+").");
		}
		//Missing/extra reps: every RECONSTRUCTED group must also appear in the sidecar.
		for(String rep : groups.keySet()){
			if(!sidecarReps.contains(rep)){
				failures.add("Group '"+rep+"' exists in cluster.tsv but has NO corresponding sidecar row (missing rep).");
			}
		}
		if(sidecarReps.size()!=groups.size()){
			failures.add("sidecar has "+sidecarReps.size()+" distinct rep rows but cluster.tsv "
				+"reconstructs "+groups.size()+" distinct groups (extra/missing rep(s)).");
		}

		result.pass=failures.isEmpty();
		return result;
	}

	/** Raw byte-stream SHA-256 of the file's literal bytes -- no line parsing, no reconstruction,
	 * so it sees exactly what any corruption (including a stray blank/partial line) would change. */
	static String hashFileRaw(String path){
		final MessageDigest md=sha256();
		try(BufferedInputStream in=new BufferedInputStream(new FileInputStream(path))){
			final byte[] buf=new byte[1<<16];
			int n;
			while((n=in.read(buf))>=0){md.update(buf, 0, n);}
		}catch(IOException e){
			throw new RuntimeException("I/O error reading "+path+" for raw hashing", e);
		}
		return toHex(md.digest());
	}

	static MessageDigest sha256(){
		try{return MessageDigest.getInstance("SHA-256");}
		catch(NoSuchAlgorithmException e){throw new RuntimeException(e);}
	}

	static String toHex(byte[] b){
		final StringBuilder sb=new StringBuilder(b.length*2);
		for(byte x : b){sb.append(String.format("%02x", x));}
		return sb.toString();
	}
}
