package prot;

import java.util.List;
import java.util.regex.Pattern;

import structures.ByteBuilder;

/**
 * Shared, independently-testable constructions used across every Candidate-C pipeline phase
 * (`CANDIDATE_C_BUILD_PIPELINE_DESIGN_v1.md` v4): the ASCII-safe key derivations, the
 * ORDER-SENSITIVE list hash (deliberately distinct from
 * {@link ArtifactRoleManifestGenerator#canonicalListSha256}'s lexically-sorted, order-BLIND
 * construction), and the HMM DATE-line canonicalization text transform. One implementation, reused
 * by every phase that needs it, rather than four subtly-different copies -- these are pure
 * UTILITY constructions with a single frozen definition each, not domain logic to independently
 * re-derive per tool (that discipline applies to things like group/member expansion, not to a
 * fixed hash recipe every consumer must compute identically by construction).
 *
 * @author Eru
 */
final class CandidateCUtil {

	private CandidateCUtil(){}

	static final Pattern ARTIFACT_KEY_PATTERN=Pattern.compile("^p_[0-9a-f]{16}$");
	static final Pattern CHALLENGE_KEY_PATTERN=Pattern.compile("^db_p_[0-9a-f]{16}$");

	/** {@code "p_" + hex16(FNV-1a-64(artifactId))} -- ASCII-safe, filesystem/HMM-NAME-safe by
	 * construction, regardless of what characters the underlying {@code family_id}/{@code rep_id}
	 * contains (this project has already seen pipe/quote-bearing rep_ids break naive shell/filename
	 * handling once). */
	static String artifactKey(final String artifactId){
		return "p_"+ArtifactRoleManifestGenerator.hex16(ReducedAlphabetSeedAssay.stableHash(artifactId));
	}

	/** {@code "db_" + <focal artifact's artifactKey>} -- the challenge-set directory/job-name key,
	 * reusing an already-safe key rather than deriving a second hash. */
	static String challengeKey(final String focalArtifactKey){
		return "db_"+focalArtifactKey;
	}

	/**
	 * The ORDER-SENSITIVE list hash used for both {@code selected_member_sha256_ordered} (§3c,
	 * over {@code selection_rank} order) and {@code expected_profile_order_sha256} (§3b, over the
	 * pinned focal-then-neighbor-rank concat order): join the IDs in the EXACT order given (never
	 * sorted), exactly one LF after every id including the last, streaming SHA-256 over that UTF-8
	 * byte sequence. Non-ASCII IDs are rejected outright (matching this project's established
	 * {@code canonicalListSha256} ASCII discipline).
	 * @param orderedIds IDs in the exact order to bind (caller-determined -- rank order, concat
	 *        order, etc; this method never reorders anything).
	 * @return Lowercase hex SHA-256.
	 */
	static String orderedIdListSha256(final List<String> orderedIds){
		final ByteBuilder bb=new ByteBuilder();
		for(final String id : orderedIds){
			for(int i=0; i<id.length(); i++){
				if(id.charAt(i)>127){
					throw new RuntimeException("ID '"+id+"' contains a non-ASCII character -- ordered "
						+"hashing requires pure-ASCII IDs.");
				}
			}
			bb.append(id).append('\n');
		}
		return ArtifactRoleManifestGenerator.sha256Hex(bb.toBytes());
	}

	/**
	 * Replaces every {@code DATE  <value>} line in a raw HMMER ASCII {@code .hmm} file's text with
	 * a fixed constant, leaving every other line (including the RNG-derived but seed-deterministic
	 * {@code STATS} block) untouched. Confirmed safe by direct measurement
	 * (`scripts/gates/candidatec_hmmer_preflight_v1.sh`, 2026-09-01): {@code hmmpress}/{@code
	 * hmmscan} accept a DATE-edited file without error, and reported bit scores are identical to
	 * the uncanonicalized file's.
	 * @param hmmText The raw {@code .hmm} file's full text (one or more concatenated profiles).
	 * @return The same text with every {@code DATE} line's value replaced.
	 */
	static String canonicalizeHmmDate(final String hmmText){
		final StringBuilder out=new StringBuilder(hmmText.length());
		final String[] lines=hmmText.split("\n", -1);
		for(int i=0; i<lines.length; i++){
			final String line=lines[i];
			if(line.startsWith("DATE  ")){out.append("DATE  [canonicalized]");}
			else{out.append(line);}
			if(i<lines.length-1){out.append('\n');}
		}
		return out.toString();
	}
}
