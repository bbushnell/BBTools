package prot;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import fileIO.ByteFile;

/**
 * Pass 0 (Increment 3B design v1-v7, slice-1 v2 correction): loads the tracked-family
 * roster canonically -- the rank&lt;-&gt;rep_id lookup Pass 0.5 and Pass 1 both need.
 * The first physical line is either the historical literal header
 * {@code #rank<TAB>rep_id<TAB>occ_total} or its post-split bound form with one final
 * {@code <TAB>roster_sha80=HEX80} token; every subsequent physical line is exactly
 * three columns {@code rank<TAB>rep_id<TAB>occ_total} -- no extra/missing columns, no
 * skippable blank or comment lines beyond that one header. {@code rep_id} is decoded with
 * strict UTF-8 (rejects invalid byte sequences, never silently substitutes); {@code rank}/
 * {@code occ_total} are canonical nonnegative decimal integers (no sign, no leading zero
 * except "0" itself). Ranks are required to be dense and contiguous from 0 -- this matches
 * the accepted family-rank model project-wide, but a real-file run still VERIFIES it here
 * rather than assuming it. {@code occ_total} is retained per rank: real roster provenance
 * and a later scale cross-check (root_review_slice1_v1.md sec1).
 *
 * <p>Small (thousands of rows), loaded once per run -- a boxed {@code HashMap<String,Integer>}
 * is appropriate here (not the true-scale hot path, which is Pass 1's member index).</p>
 *
 * @author Eru
 */
public final class HbmRosterLoader {

	private static final String HEADER_LINE="#rank\trep_id\tocc_total";
	private static final String BOUND_HEADER_PREFIX=HEADER_LINE+"\troster_sha80=";

	public final String[] repIdByRank;
	public final long[] occTotalByRank;
	/** Accepted post-split roster binding, or null for the historical unbound schema. */
	public final String sourceRosterSha80;
	private final String headerLine;
	private final HashMap<String,Integer> rankByRepId;

	private HbmRosterLoader(final String[] repIdByRank, final long[] occTotalByRank, final HashMap<String,Integer> rankByRepId,
			final String sourceRosterSha80_, final String headerLine_){
		this.repIdByRank=repIdByRank;
		this.occTotalByRank=occTotalByRank;
		this.rankByRepId=rankByRepId;
		sourceRosterSha80=sourceRosterSha80_; headerLine=headerLine_;
	}

	public int size(){return repIdByRank.length;}

	/** Returns the rank for repId, or -1 if repId is not a tracked roster representative. */
	public int rankOf(final String repId){
		final Integer r=rankByRepId.get(repId);
		return r==null ? -1 : r.intValue();
	}

	public boolean legalRank(final int rank){return rank>=0 && rank<repIdByRank.length;}

	public static HbmRosterLoader load(final String fname){
		final ByteFile bf=ByteFile.makeByteFile(fname, true);
		try{
			final ArrayList<String> reps=new ArrayList<String>();
			final ArrayList<Long> occTotals=new ArrayList<Long>();
			final HashMap<String,Integer> rankByRepId=new HashMap<String,Integer>();
			long lineNum=0;
			int expectedRank=0;
			boolean sawHeader=false; String sourceRosterSha80=null,headerLine=null;
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lineNum++;
				if(!sawHeader){
					final String headerText=HbmCanonicalParse.utf8StrictDecode(line);
					if(HEADER_LINE.equals(headerText)){headerLine=HEADER_LINE;}
					else if(headerText.startsWith(BOUND_HEADER_PREFIX)){
						sourceRosterSha80=headerText.substring(BOUND_HEADER_PREFIX.length());
						HbmCanonicalParse.checkLowercaseHexSha80(sourceRosterSha80,"execution roster source roster sha80");
						headerLine=headerText;
					}else{throw new RuntimeException("ROSTER_HEADER: roster line 1 must be the historical header or its roster-bound extension, got '"+headerText+"'");}
					sawHeader=true;
					continue;
				}
				final int tab1=HbmCanonicalParse.indexOf(line, (byte)'\t', 0);
				if(tab1<0){throw new RuntimeException("ROSTER_SHAPE: roster row "+lineNum+" is missing columns (no TAB): '"+HbmCanonicalParse.utf8StrictDecode(line)+"'");}
				final int tab2=HbmCanonicalParse.indexOf(line, (byte)'\t', tab1+1);
				if(tab2<0){throw new RuntimeException("ROSTER_SHAPE: roster row "+lineNum+" has fewer than 3 columns: '"+HbmCanonicalParse.utf8StrictDecode(line)+"'");}
				final int tab3=HbmCanonicalParse.indexOf(line, (byte)'\t', tab2+1);
				if(tab3>=0){throw new RuntimeException("ROSTER_SHAPE: roster row "+lineNum+" has more than 3 columns: '"+HbmCanonicalParse.utf8StrictDecode(line)+"'");}

				final String rankStr=HbmCanonicalParse.utf8StrictDecode(line, 0, tab1);
				final String repId=HbmCanonicalParse.utf8StrictDecode(line, tab1+1, tab2-tab1-1);
				final String occStr=HbmCanonicalParse.utf8StrictDecode(line, tab2+1, line.length-tab2-1);

				final long rank=HbmCanonicalParse.parseCanonicalNonNegLong(rankStr, "roster rank at row "+lineNum);
				final long occTotal=HbmCanonicalParse.parseCanonicalNonNegLong(occStr, "roster occ_total at row "+lineNum);
				if(rank!=expectedRank){
					throw new RuntimeException("ROSTER_RANK_DENSE: roster rank not dense/contiguous from 0 at row "+lineNum+
						": expected "+expectedRank+" got "+rank);
				}
				if(repId.isEmpty()){throw new RuntimeException("ROSTER_SHAPE: roster row "+lineNum+" (rank "+rank+") has an empty rep_id");}
				if(rankByRepId.put(repId, Integer.valueOf((int)rank))!=null){
					throw new RuntimeException("ROSTER_DUP_REPID: duplicate roster rep_id '"+repId+"' at row "+lineNum);
				}
				reps.add(repId);
				occTotals.add(Long.valueOf(occTotal));
				expectedRank++;
			}
			if(!sawHeader){throw new RuntimeException("ROSTER_HEADER: roster file has no header line: "+fname);}
			if(reps.isEmpty()){throw new RuntimeException("ROSTER_EMPTY: roster has zero data rows: "+fname);}
			final long[] occArr=new long[occTotals.size()];
			for(int i=0; i<occArr.length; i++){occArr[i]=occTotals.get(i).longValue();}
			final HbmRosterLoader result=new HbmRosterLoader(reps.toArray(new String[reps.size()]), occArr, rankByRepId,
				sourceRosterSha80,headerLine);
			requireByteCanonical(fname, result);
			return result;
		}finally{
			bf.close();
		}
	}

	/** {@link ByteFile}'s line splitting silently strips a trailing CR immediately before an
	 *  LF when it finds one (slice-1 v6, root_review_slice1_v5.md blocker 1: root demonstrated
	 *  a well-formed CRLF roster parses successfully through the loop above, identically to its
	 *  LF-only equivalent) -- so nothing in the per-line parsing above can tell a CRLF file from
	 *  a bare-LF file. That means {@link #toCanonicalBytes()}'s LF-only reconstruction is not
	 *  provably the hash of the bytes {@link #load} actually accepted for such a file. Close
	 *  this the way root's own review suggests as simplest (the roster is small, per the class
	 *  doc above): re-read the exact decoded file bytes, bypassing line normalization, and
	 *  require them to already be byte-for-byte identical to {@link #toCanonicalBytes()}'s
	 *  reconstruction -- rejecting CRLF, a missing terminal LF, or any other non-canonical
	 *  physical representation outright, rather than trying to reproduce whatever convention was
	 *  present. */
	private static void requireByteCanonical(final String fname, final HbmRosterLoader result){
		final byte[] raw;
		try{raw=MagQCTextResource.bytes(fname);}
		catch(IOException e){throw new RuntimeException("ROSTER_IO: failed to re-read raw bytes for canonical verification: "+fname, e);}
		final byte[] canonical=result.toCanonicalBytes();
		if(Arrays.equals(raw, canonical)){return;}
		for(final byte b : raw){
			if(b=='\r'){throw new RuntimeException("ROSTER_CR_REJECTED: roster file contains a carriage return byte; only bare-LF line endings are accepted: "+fname);}
		}
		if(raw.length==0 || raw[raw.length-1]!='\n'){
			throw new RuntimeException("ROSTER_MISSING_TERMINAL_LF: roster file does not end with a single LF byte: "+fname);
		}
		throw new RuntimeException("ROSTER_NONCANONICAL_BYTES: roster file's raw bytes are not byte-equal to its canonical LF-only serialization: "+fname);
	}

	/** Reconstructs the exact canonical roster-file bytes this object would have been parsed
	 *  from -- valid because {@link #load} enforces a strict canonical grammar (exact header,
	 *  dense ranks from 0, canonical decimal fields, no blank/comment lines), so there is only
	 *  one canonical byte serialization for a given (rank,rep_id,occ_total) sequence. */
	public byte[] toCanonicalBytes(){
		final StringBuilder sb=new StringBuilder(headerLine).append('\n');
		for(int i=0; i<repIdByRank.length; i++){
			sb.append(i).append('\t').append(repIdByRank[i]).append('\t').append(occTotalByRank[i]).append('\n');
		}
		return sb.toString().getBytes(StandardCharsets.UTF_8);
	}

	/** Requires this roster object's canonical reserialization to hash to a TRUSTED, previously
	 *  and independently verified {@code roster_whole_file_sha256} value (slice-1 v3,
	 *  root_review_slice1_v3.md sec1/2: verifies the OBJECT's bytes, not merely a caller-supplied
	 *  path string, against a value the caller already trusts). This is what catches a
	 *  same-size roster object whose rank<->rep_id mapping has been reordered or altered --
	 *  reserializing it produces different bytes and therefore a different hash. */
	public void requireMatchesTrustedHash(final String trustedRosterWholeFileSha256, final String tag){
		final byte[] bytes=toCanonicalBytes();
		final String recomputed=HbmMemberIndexFormat.toHexLower(HbmMemberIndexFormat.sha256(bytes, 0, bytes.length));
		if(!recomputed.equals(trustedRosterWholeFileSha256)){
			throw new RuntimeException(tag+": roster object's canonical reserialization hash "+recomputed+
				" does not match the trusted roster_whole_file_sha256 "+trustedRosterWholeFileSha256);
		}
	}

	/** Sha80 textual-provenance counterpart to {@link #requireMatchesTrustedHash}. The
	 *  roster's canonical bytes are still hashed with SHA-256; only the final 80 bits are
	 *  represented in the new text contract. */
	public void requireMatchesTrustedSha80(final String trustedRosterWholeFileSha80, final String tag){
		HbmCanonicalParse.checkLowercaseHexSha80(trustedRosterWholeFileSha80, "trusted roster sha80");
		final String recomputed=DigestSuffix.bytes(toCanonicalBytes());
		if(!recomputed.equals(trustedRosterWholeFileSha80)){
			throw new RuntimeException(tag+": roster object's canonical reserialization sha80 "+recomputed+
				" does not match the trusted roster_whole_file_sha80 "+trustedRosterWholeFileSha80);
		}
	}

	/** Requires two roster objects to describe exactly the same rank/rep_id/occ_total content
	 *  (slice-1 v2/v3: "compare every rank/rep/occurrence byte-equivalently" -- used wherever a
	 *  caller-supplied roster object must be proven equivalent to one derived from a
	 *  hash-checked file, rather than trusted as-is). */
	public static void requireEquivalent(final HbmRosterLoader a, final HbmRosterLoader b, final String tag){
		if(a.size()!=b.size()){
			throw new RuntimeException(tag+": roster objects have different sizes ("+a.size()+" vs "+b.size()+")");
		}
		for(int i=0; i<a.size(); i++){
			if(!a.repIdByRank[i].equals(b.repIdByRank[i])){
				throw new RuntimeException(tag+": roster objects differ at rank "+i+"'s rep_id");
			}
			if(a.occTotalByRank[i]!=b.occTotalByRank[i]){
				throw new RuntimeException(tag+": roster objects differ at rank "+i+"'s occ_total");
			}
		}
	}
}
