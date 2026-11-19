package prot;

import java.util.Arrays;

import ddl.DDLFamilyProfiles;
import ddl.DDLKeyTieBreak;
import ddl.DDLPentamerBlacklist;
import ddl.ProteinKmerGeometry;
import ddl.ProteinPositionMode;
import structures.ByteBuilder;

/** Prepares one family's member-derived payload without intermediate files or negative calibration.
 * A worker is single-threaded; parallel callers each own one worker. Members must be supplied in
 * canonical order and remain unchanged during the call. Results retain no member arrays or worker scratch.
 * @author Keqing
 */
public final class FamilyMemberPreparation{

	public FamilyMemberPreparation(int partitions_, int overlap_, DDLPentamerBlacklist blacklist_,
			FamilyCalibration.MemberRetention retention_){
		this(partitions_, overlap_, blacklist_, retention_, false);
	}

	public FamilyMemberPreparation(int partitions_, int overlap_, DDLPentamerBlacklist blacklist_,
			FamilyCalibration.MemberRetention retention_, boolean retainProfileRows_){
		this(ProteinKmerGeometry.exact(5), partitions_, overlap_, blacklist_, retention_, retainProfileRows_);
	}

	public FamilyMemberPreparation(ProteinKmerGeometry geometry_, int partitions_, int overlap_,
			DDLPentamerBlacklist blacklist_, FamilyCalibration.MemberRetention retention_, boolean retainProfileRows_){
		this(geometry_, partitions_, overlap_, blacklist_, retention_, retainProfileRows_, MemberLengthRange.LEGACY);
	}

	public FamilyMemberPreparation(ProteinKmerGeometry geometry_, int partitions_, int overlap_,
			DDLPentamerBlacklist blacklist_, FamilyCalibration.MemberRetention retention_, boolean retainProfileRows_, MemberLengthRange memberRange_){
		this(geometry_, partitions_, overlap_, blacklist_, retention_, retainProfileRows_, memberRange_, ProteinPositionMode.FIXED);
	}

	/** One coordinate rule governs every retained member; the plain reference row is also retained. */
	public FamilyMemberPreparation(ProteinKmerGeometry geometry_, int partitions_, int overlap_,
			DDLPentamerBlacklist blacklist_, FamilyCalibration.MemberRetention retention_, boolean retainProfileRows_,
			MemberLengthRange memberRange_, ProteinPositionMode positionMode_){
		this(geometry_, partitions_, overlap_, blacklist_, retention_, retainProfileRows_, memberRange_, positionMode_, FamilyAlignmentMetric.LEGACY);
	}

	/** Explicit metric selection requires matching negative calibration and export before publication. */
	public FamilyMemberPreparation(ProteinKmerGeometry geometry_, int partitions_, int overlap_,
			DDLPentamerBlacklist blacklist_, FamilyCalibration.MemberRetention retention_, boolean retainProfileRows_,
			MemberLengthRange memberRange_, ProteinPositionMode positionMode_, FamilyAlignmentMetric alignmentMetric_){
		this(geometry_, partitions_, overlap_, blacklist_, retention_, retainProfileRows_, memberRange_, positionMode_, alignmentMetric_, 1024);
	}

	/** Explicit width applies before member voting; folding an already merged profile would change its modal keys. */
	public FamilyMemberPreparation(ProteinKmerGeometry geometry_, int partitions_, int overlap_,
			DDLPentamerBlacklist blacklist_, FamilyCalibration.MemberRetention retention_, boolean retainProfileRows_,
			MemberLengthRange memberRange_, ProteinPositionMode positionMode_, FamilyAlignmentMetric alignmentMetric_, int width_){
		this(geometry_,partitions_,overlap_,blacklist_,retention_,retainProfileRows_,memberRange_,positionMode_,alignmentMetric_,width_,null);
	}

	/** Only the selected positional row uses this context; scalar/alignment measurements are unchanged. */
	public FamilyMemberPreparation(ProteinKmerGeometry geometry_,int partitions_,int overlap_,
			DDLPentamerBlacklist blacklist_,FamilyCalibration.MemberRetention retention_,boolean retainProfileRows_,
			MemberLengthRange memberRange_,ProteinPositionMode positionMode_,FamilyAlignmentMetric alignmentMetric_,int width_,
			DDLKeyTieBreak.PolicyContext keyPolicy_){
		require(width_==512 || width_==1024, "Member preparation supports only the validated 512/1024-register kernels");
		width=width_; sketch=new char[width];
		require(alignmentMetric_!=null, "Member preparation must name the alignment contract used by its lower gates");
		alignmentMetric=alignmentMetric_;
		alignmentWorker=alignmentMetric.same(FamilyAlignmentMetric.LEGACY) ? null : alignmentMetric.newWorker();
		require(memberRange_!=null, "Member preparation requires an explicit length contract");
		require(positionMode_!=null, "Member preparation requires an explicit position mode");
		memberRange=memberRange_;
		positionMode=positionMode_;
		if(geometry_==null){throw new IllegalArgumentException("Member preparation requires a protein k-mer geometry");}
		if(positionMode.positional()){DDLFamilyProfiles.validatePositionalLayout(partitions_, overlap_, geometry_,positionMode_);}
		require(retention_!=null, "Member preparation requires an explicit lower-retention policy");
		geometry=geometry_; partitions=partitions_; overlap=overlap_; blacklist=ddl.ProteinBlacklistPolicy.resolve(geometry_,blacklist_); retention=retention_;
		FamilyBlacklistBinding.validateGeometry(blacklist,geometry);
		blacklistSha80=blacklist==null ? "none" : blacklist.sha80();
		retainProfileRows=retainProfileRows_;
		keyPolicy=keyPolicy_;
		if(keyPolicy!=null){
			keyPolicy.requireGeometry(width,geometry,0,positionMode,partitions,overlap,blacklist==null ? null : blacklist.sha80());
			require(positionMode.positional() || keyPolicy.policy==DDLKeyTieBreak.Policy.SMALLEST && keyPolicy.consensusMemberBonus==0,
				"Nondefault preparation and positive consensus bonus require the selected positional row");
		}
		keyScratch=keyPolicy==null ? null : new DDLKeyTieBreak.Workspace(width);
		consensusKeys=keyPolicy!=null && (keyPolicy.policy==DDLKeyTieBreak.Policy.CONSENSUS_DIVERSITY || keyPolicy.consensusMemberBonus>0) ?
			new DDLKeyTieBreak.KeySet(width) : null;
	}

	/** Two in-memory passes preserve scalar order and score every member against the complete dimer aggregate.
	 * Input tails obey the configured member range (default100..200); the consensus is bounded by1000 residues.
	 * Empty families have no fitted payload and are explicit through members==0 and null payload fields.
	 */
	public Prepared prepare(int rank, String id, byte[][] members, byte[] consensus){
		require(rank>=0 && id!=null && id.matches("[A-Za-z0-9_.-]+") && members!=null,
			"Member preparation requires a canonical family identity and member array");
		require(consensus!=null && consensus.length>0 && consensus.length<=1000, "Fixed consensus must contain1..1000 residues");
		final int n=members.length;
		if(n==0){
			return new Prepared(rank, id, 0, 0, 0, null, null, null, null, null, null, null, null, null, geometry.name(), partitions, overlap, blacklistSha80, positionMode, alignmentMetric, null, null, width, blacklist);
		}
		final ProteinSequence reference=new ProteinSequence(id, consensus);
		require(reference.length()==consensus.length, "Prepared consensus must not contain a terminal stop marker");
		final DDLFamilyProfiles.Profile votes=new DDLFamilyProfiles.Profile(width);
		final ProteinSpectrum.Spectrum spectrum=new ProteinSpectrum.Spectrum(2);
		final FamilyMemberScalars.Moments moments=new FamilyMemberScalars.Moments(1);
		long residues=0, words=0;
		for(byte[] member:members){
			if(member==null || !memberRange.contains(member.length)){
				throw new IllegalArgumentException("Original member length must be in "+memberRange+": "+id);
			}
			final int valid=ConsensusScalars.computeScalars(member, member.length, scalars, 0, scalarMask);
			require(valid==member.length && scalarMask[0]==ConsensusScalars.ALL_DEFINED,
				"Member statistics require canonical amino acids and all five measured scalars: "+id);
			moments.add(0, scalars); spectrum.add(member);
			residues=Math.addExact(residues, member.length);
			words=Math.addExact(words, DDLFamilyProfiles.sketch(member, geometry, 0, sketch, blacklist)); votes.add(sketch);
		}
		spectrum.validate();
		require(spectrum.records==n && moments.n[0]==n && spectrum.residues==residues, "Every member must contribute to the same family population");
		final ProfileRow plain=modal(votes, n, false);
		votes.clear();
		final double norm=ProteinSpectrum.norm(spectrum);
		final double[][] scores=new double[3][n];
		final int[] alignmentCounts=new int[ALIGNMENT_STATES];
		int measured=0;
		long positionalWords=0;
		for(int i=0; i<n; i++){
			final byte[] member=members[i];
			positionalWords=Math.addExact(positionalWords, positionMode.positional() ?
				DDLFamilyProfiles.sketchPositional(member, geometry, 0, partitions, overlap, positionMode, sketch, blacklist) :
				DDLFamilyProfiles.sketch(member, geometry, 0, sketch, blacklist));
			votes.add(sketch);
			querySpectrum.clear(); querySpectrum.add(member);
			final ProteinSequence query=buffers.encode(member);
			if(alignmentWorker==null){
				// Keep the historical arithmetic and full-member score arrays byte-compatible.
				FamilyCalibration.measure(query, querySpectrum, reference, spectrum, norm, metrics);
				scores[0][i]=metrics.cosine; scores[1][measured]=metrics.identityPct; scores[2][measured++]=metrics.blosum;
				alignmentCounts[FamilyAlignmentMetric.State.MEASURED.ordinal()]++;
			}else{
				final double cosine=ProteinSpectrum.cosine(querySpectrum, spectrum, ProteinSpectrum.norm(querySpectrum), norm);
				require(Double.isFinite(cosine) && cosine>=0 && cosine<=1, "Member cosine remains a measured full-sequence metric");
				scores[0][i]=cosine;
				alignmentWorker.measure(query, reference, ProteinSearcher.selfScore(query), alignment);
				alignmentCounts[alignment.state.ordinal()]++;
				if(alignment.measured()){
					scores[1][measured]=alignment.identityPct(); scores[2][measured++]=alignment.blosum();
				}
			}
		}
		require(positionalWords==words, "Plain and positional profiles must count the same original valid K5 windows");
		final ProfileRow positional=modal(votes, n, retainProfileRows);
		if(keyPolicy!=null){
			if(consensusKeys!=null){
				DDLFamilyProfiles.enumerateKeysEncoded(reference.enc,geometry,0,positionMode,partitions,overlap,blacklist,consensusKeys);
			}
			final boolean hash=keyPolicy.policy==DDLKeyTieBreak.Policy.DIVERSITY || keyPolicy.policy==DDLKeyTieBreak.Policy.CONSENSUS_DIVERSITY;
			votes.writePackedValues(keyPolicy,hash ? DDLKeyTieBreak.familySeed(id) : 0,consensusKeys,keyScratch,null,positional.values);
			if(keyPolicy.consensusMemberBonus>0 && positional.support!=null){
				for(int bucket=0;bucket<width;bucket++){
					final int key=positional.values[bucket];
					positional.support[bucket]=key==0 ? 0 : votes.memberSupport(bucket,key);
				}
			}
		}
		final FamilyCalibration.Quantile[] lower=alignmentWorker==null ? FamilyCalibration.memberCutoffs(scores, retention, null) : null;
		final FamilyCalibration.MemberCutoff[] measuredLower=alignmentWorker==null ? null :
			FamilyCalibration.measuredMemberCutoffs(scores, new int[]{n, measured, measured}, n, retention);
		final double[] statistics=new double[20]; moments.exportStatistics(0, statistics, decimal);
		return new Prepared(rank, id, n, residues, words, Arrays.copyOf(consensus, consensus.length), reference,
			plain.values, positional.values, positional.support, positional.occupancy, spectrum, lower, statistics, geometry.name(), partitions, overlap, blacklistSha80, positionMode,
			alignmentMetric, measuredLower, alignmentCounts, width, blacklist);
	}

	private ProfileRow modal(DDLFamilyProfiles.Profile votes, int members, boolean retain){
		require(votes.members()==members && votes.buckets()==width,
			"Modal votes must include every original member exactly once at the worker's declared width");
		final char[] values=new char[width];
		final int[] support=retain ? new int[width] : null, occupancy=retain ? new int[width] : null;
		for(int i=0; i<values.length; i++){
			values[i]=votes.value(i, false);
			if(retain){support[i]=votes.votes(i, false); occupancy[i]=votes.occupied(i);}
		}
		return new ProfileRow(values, support, occupancy);
	}

	/** Caller-owned detached payload. Treat its arrays as read-only while sharing with downstream workers.
	 * Raw long dimer counts are retained for the legacy negative-cosine maximum; float projection comes later.
	 */
	public static final class Prepared{
		private Prepared(int rank_, String id_, int members_, long residues_, long words_, byte[] raw_, ProteinSequence consensus_,
				char[] plain_, char[] positional_, int[] positionalSupport_, int[] positionalOccupancy_, ProteinSpectrum.Spectrum spectrum_,
				FamilyCalibration.Quantile[] lower_, double[] statistics_, String geometry_, int partitions_, int overlap_, String blacklistSha80_,
				ProteinPositionMode positionMode_, FamilyAlignmentMetric alignmentMetric_,
				FamilyCalibration.MemberCutoff[] measuredLower_, int[] alignmentCounts_, int width_, DDLPentamerBlacklist blacklist_){
			assert(width_==512 || width_==1024) : "Prepared profile width must name a validated member-sketch kernel";
			assert(members_==0 ? plain_==null && positional_==null :
				plain_!=null && positional_!=null && plain_.length==width_ && positional_.length==width_) :
				"Downstream calibration requires both modal profiles at the declared width, or neither for an empty family";
			assert(alignmentMetric_!=null) : "Prepared gates must retain their exact metric descriptor";
			assert(members_==0 ? spectrum_==null && lower_==null && measuredLower_==null && statistics_==null :
				spectrum_!=null && (alignmentMetric_.same(FamilyAlignmentMetric.LEGACY) ?
					lower_!=null && lower_.length==3 && measuredLower_==null :
					lower_==null && measuredLower_!=null && measuredLower_.length==3) && statistics_!=null && statistics_.length==20) :
				"Downstream calibration needs a complete fitted payload or an explicit empty family";
			if(members_>0){
				require(alignmentCounts_!=null && alignmentCounts_.length==ALIGNMENT_STATES,
					"Each original member needs exactly one explicit alignment outcome");
				int total=0; for(int count:alignmentCounts_){require(count>=0, "Alignment outcome counts cannot be negative"); total=Math.addExact(total,count);}
				require(total==members_, "Alignment outcomes must conserve the full member population");
			}
			assert((positionalSupport_==null)==(positionalOccupancy_==null)) :
				"Positional DDL profile rows require real support and occupancy arrays or neither";
			assert(positionalSupport_==null || positionalSupport_.length==width_ && positionalOccupancy_.length==width_) :
				"Serialized modal support and occupancy must cover exactly the declared profile width";
			rank=rank_; id=id_; members=members_; residues=residues_; validWords=words_; consensusRaw=raw_; consensus=consensus_;
			plain=plain_; positional=positional_; positionalSupport=positionalSupport_; positionalOccupancy=positionalOccupancy_;
			spectrum=spectrum_; lower=lower_; statistics=statistics_;
			geometry=geometry_; partitions=partitions_; overlap=overlap_; blacklistSha80=blacklistSha80_;
			blacklist=blacklist_;
			require(blacklistSha80.equals(blacklist==null ? "none" : blacklist.sha80()),
				"Prepared profile identity must describe the filter instance used for member sketches");
			positionMode=positionMode_;
			width=width_;
			alignmentMetric=alignmentMetric_; measuredLower=measuredLower_; alignmentCounts=alignmentCounts_;
		}

		/** Appends the selected replacement-profile row without clearing; false for empty families.
		 * The historical method name also covers the explicit position-free control.
		 */
		public boolean appendPositionalProfile(ByteBuilder line){
			if(members==0){return false;}
			require(positional!=null && positionalSupport!=null && positionalOccupancy!=null,
				"Positional profile output requires opt-in support and occupancy retention");
			line.append("tail\t").append(id).tab().append(width).tab().append(members).tab().append(residues).tab().append(validWords).append("\tnonempty");
			append(line, positional); append(line, positionalSupport); append(line, positionalOccupancy);
			return true;
		}

		/** A direct resident handoff shares the object; only measured legacy reloads may use an equivalent instance. */
		public void verifyBlacklist(DDLPentamerBlacklist expected,boolean direct){
			if(blacklist==expected){return;}
			require(!direct && blacklist!=null && expected!=null && !expected.sha80().equals("unmeasured") &&
				blacklistSha80.equals(expected.sha80()) && blacklist.wordLength()==expected.wordLength() &&
				blacklist.matchMode()==expected.matchMode() && blacklist.size()==expected.size() &&
				blacklist.provenance().equals(expected.provenance()),
				"Prepared members and resident queries must share the same unmeasured filter object: "+id);
		}

		public final int rank, members;
		public final String id;
		/** Both modal rows share this geometry and blacklist; only the positional row uses the layout. */
		public final String geometry, blacklistSha80;
		/** Shared immutable filter; an unmeasured label cannot identify its word set. */
		public final DDLPentamerBlacklist blacklist;
		public final ProteinPositionMode positionMode;
		public final int partitions, overlap, width;
		public final long residues, validWords;
		public final byte[] consensusRaw;
		public final ProteinSequence consensus;
		/** positional carries the selected replacement layout identified by positionMode. */
		public final char[] plain, positional;
		public final int[] positionalSupport, positionalOccupancy;
		public final ProteinSpectrum.Spectrum spectrum;
		public final FamilyCalibration.Quantile[] lower;
		/** Legacy cutoffs above are absent for explicit metrics; these retain assigned/failed counts and unavailable state. */
		public final FamilyCalibration.MemberCutoff[] measuredLower;
		public final FamilyAlignmentMetric alignmentMetric;
		/** Counts by FamilyAlignmentMetric.State ordinal; null only for a zero-member family. */
		public final int[] alignmentCounts;
		public final double[] statistics;
	}

	private static void append(ByteBuilder line, char[] values){
		line.tab(); for(int i=0; i<values.length; i++){if(i>0){line.append(',');} line.append((int)values[i]);}
	}

	private static void append(ByteBuilder line, int[] values){
		line.tab(); for(int i=0; i<values.length; i++){if(i>0){line.append(',');} line.append(values[i]);}
	}

	private static final class ProfileRow{
		private ProfileRow(char[] values_, int[] support_, int[] occupancy_){
			values=values_; support=support_; occupancy=occupancy_;
		}
		final char[] values;
		final int[] support, occupancy;
	}

	private static void require(boolean condition, String reason){if(!condition){throw new IllegalArgumentException(reason);}}
	private final DDLKeyTieBreak.PolicyContext keyPolicy;
	private final DDLKeyTieBreak.Workspace keyScratch;
	private final DDLKeyTieBreak.KeySet consensusKeys;
	private final ProteinKmerGeometry geometry;
	private final ProteinPositionMode positionMode;
	private final MemberLengthRange memberRange;
	private final int partitions, overlap, width;
	private final DDLPentamerBlacklist blacklist;
	private final String blacklistSha80;
	private final FamilyCalibration.MemberRetention retention;
	private final FamilyAlignmentMetric alignmentMetric;
	private final FamilyAlignmentMetric.Worker alignmentWorker;
	private final FamilyAlignmentMetric.Observation alignment=new FamilyAlignmentMetric.Observation();
	private final boolean retainProfileRows;
	private final FamilyCalibration.Buffers buffers=new FamilyCalibration.Buffers();
	private final ProteinSpectrum.Spectrum querySpectrum=new ProteinSpectrum.Spectrum(2);
	private final FamilyCalibration.Metrics metrics=new FamilyCalibration.Metrics();
	private final char[] sketch;
	private final double[] scalars=new double[5];
	private final int[] scalarMask=new int[1];
	private final ByteBuilder decimal=new ByteBuilder();
	private static final int ALIGNMENT_STATES=FamilyAlignmentMetric.State.values().length;
}
