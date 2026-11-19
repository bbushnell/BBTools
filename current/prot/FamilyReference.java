package prot;

import dna.AminoAcid;

/** Published family data. Arrays are transferred to this object and must not be mutated.
 * @author Sayu
 */
public final class FamilyReference{

	/** Empty roster entries have no gates or payloads; their unused numeric fields are finite zeros. */
	public static FamilyReference empty(String id,int rank){
		return new FamilyReference(id,rank,0,0,0,null,null,null,0,0,0,0,0,0);
	}

	public FamilyReference(String id,int rank,long members,long residues,long validWords,
			ProteinSequence consensus,char[] ddl,float[] spectrum,float cosineLower,float identityLower,
			float blosumLower,float cosineAccept,float identityUpper,float blosumFinal){
		this(id,rank,members,residues,validWords,consensus,ddl,spectrum,cosineLower,identityLower,
			blosumLower,cosineAccept,identityUpper,blosumFinal,raw(consensus));
	}

	/** Raw bytes preserve ambiguous letters/case that ProteinSequence intentionally folds. */
	public FamilyReference(String id,int rank,long members,long residues,long validWords,
			ProteinSequence consensus,char[] ddl,float[] spectrum,float cosineLower,float identityLower,
			float blosumLower,float cosineAccept,float identityUpper,float blosumFinal,byte[] consensusRaw){
		this(id,rank,members,residues,validWords,consensus,ddl,spectrum,cosineLower,identityLower,
			blosumLower,cosineAccept,identityUpper,blosumFinal,consensusRaw,0,new double[10]);
	}

	/** Scalar bounds are interleaved inclusive minima/maxima; a clear mask bit requires two inert zeros. */
	public FamilyReference(String id,int rank,long members,long residues,long validWords,
			ProteinSequence consensus,char[] ddl,float[] spectrum,float cosineLower,float identityLower,
			float blosumLower,float cosineAccept,float identityUpper,float blosumFinal,byte[] consensusRaw,
			int scalarMask,double[] scalarBounds){
		this(id,rank,members,residues,validWords,consensus,ddl,spectrum,cosineLower,identityLower,
			blosumLower,cosineAccept,identityUpper,blosumFinal,consensusRaw,scalarMask,scalarBounds,
			AlignmentCalibration.ENABLED,false);
	}

	/** Explicitly records whether a live family has measured identity/BLOSUM member cutoffs. */
	public FamilyReference(String id,int rank,long members,long residues,long validWords,
			ProteinSequence consensus,char[] ddl,float[] spectrum,float cosineLower,float identityLower,
			float blosumLower,float cosineAccept,float identityUpper,float blosumFinal,byte[] consensusRaw,
			int scalarMask,double[] scalarBounds,AlignmentCalibration alignmentCalibration){
		this(id,rank,members,residues,validWords,consensus,ddl,spectrum,cosineLower,identityLower,
			blosumLower,cosineAccept,identityUpper,blosumFinal,consensusRaw,scalarMask,scalarBounds,alignmentCalibration,true);
	}

	private FamilyReference(String id,int rank,long members,long residues,long validWords,
			ProteinSequence consensus,char[] ddl,float[] spectrum,float cosineLower,float identityLower,
			float blosumLower,float cosineAccept,float identityUpper,float blosumFinal,byte[] consensusRaw,
			int scalarMask,double[] scalarBounds,AlignmentCalibration alignmentCalibration,
			boolean explicitAlignmentCalibration){
		require(alignmentCalibration!=null,"Missing alignment calibration availability: "+id);
		this.alignmentCalibration=alignmentCalibration; this.explicitAlignmentCalibration=explicitAlignmentCalibration;
		require(scalarMask>=0 && scalarMask<32 && scalarBounds!=null && scalarBounds.length==10,
			"Scalar gates require a five-bit mask and ten bounds");
		for(int i=0;i<5;i++){
			final double lo=scalarBounds[2*i],hi=scalarBounds[2*i+1];
			require(Double.isFinite(lo) && Double.isFinite(hi),"Scalar bounds must be finite: "+id);
			require((scalarMask&(1<<i))!=0 ? members>=30 && lo<hi : lo==0 && hi==0,
				"Enabled scalar gates require N>=30 and increasing bounds; disabled gates require zeros: "+id);
		}
		this.scalarMask=scalarMask; this.scalarBounds=scalarBounds;
		require(id!=null && !id.isEmpty(),"Missing family ID");
		for(int i=0;i<id.length();i++){
			require(id.charAt(i)>32 && id.charAt(i)<127,"Family ID must be a printable ASCII token");
		}
		require(rank>=0 && members>=0 && residues>=0 && validWords>=0 && validWords<=residues,"Invalid family metadata: "+id);
		require(Float.isFinite(cosineLower) && Float.isFinite(identityLower) && Float.isFinite(blosumLower) &&
			Float.isFinite(cosineAccept) && Float.isFinite(identityUpper) && Float.isFinite(blosumFinal),
			"Every reference numeric field must be finite, including empty families: "+id);
		this.id=id; this.rank=rank; this.members=members; this.residues=residues; this.validWords=validWords;
		this.consensus=consensus; this.consensusRaw=consensusRaw; this.ddl=ddl; this.spectrum=spectrum;
		this.cosineLower=cosineLower; this.identityLower=identityLower; this.blosumLower=blosumLower;
		this.cosineAccept=cosineAccept; this.identityUpper=identityUpper; this.blosumFinal=blosumFinal;
		if(members==0){
			require(residues==0 && validWords==0 && consensus==null && consensusRaw==null && ddl==null && spectrum==null,
				"Empty family must have zero metadata and null payloads: "+id);
			require(cosineLower==0 && identityLower==0 && blosumLower==0 && cosineAccept==0 &&
				identityUpper==0 && blosumFinal==0,
				"Empty family has no measured gates; unused numeric fields must be zero: "+id);
			require(alignmentCalibration==AlignmentCalibration.ENABLED,"Empty family has no alignment calibration channel: "+id);
			spectrumNorm=0;
		}else{
			require(consensus!=null && id.equals(consensus.id) && consensusRaw!=null && consensusRaw.length>0,
				"Missing or mismatched consensus: "+id);
			require(consensusRaw.length==consensus.enc.length,"Raw/encoded consensus length disagreement: "+id);
			for(int i=0;i<consensusRaw.length;i++){
				final int c=consensusRaw[i]&255;
				// Same accepted alphabet as Blosum62.encodeResidue; no second encoded array.
				final int standard=c<128 ? AminoAcid.acidToNumber[c] : -1,upper=c&~32;
				final int encoded=standard>=0 && standard<20 ? standard :
					upper=='X' || upper=='B' || upper=='Z' || upper=='J' ? Blosum62.X_CODE : -1;
				if(encoded<0 || encoded!=consensus.enc[i]){
					throw new IllegalArgumentException("Raw/encoded consensus disagreement: "+id+" at "+i);
				}
			}
			require(ddl!=null && FamilyReferenceCodec.validDdlWidth(ddl.length) && spectrum!=null &&
				(spectrum.length==FamilyReferenceCodec.SPECTRUM_WIDTH || spectrum.length==FamilyReferenceCodec.MONOMER_WIDTH),
				"Reference requires 512/1024 DDL registers and 400 dimers or 20 projected monomers: "+id);
			require(Float.isFinite(cosineLower) && cosineLower>=0 && cosineLower<=1 &&
				Float.isFinite(cosineAccept) && (double)cosineAccept>=0.72 && cosineAccept>=cosineLower,
				"Invalid cosine gates: "+id);
			if(alignmentCalibration==AlignmentCalibration.ENABLED){
				require(Float.isFinite(identityLower) && identityLower>=0 && identityLower<=100 && Float.isFinite(blosumLower),
					"Invalid lower alignment gates: "+id);
				require(upper(identityUpper,Math.nextUp(110f)) && Float.isFinite(blosumFinal) &&
					blosumFinal>=0 && blosumFinal>=blosumLower,"Invalid upper/final alignment gates: "+id);
				require(identityUpper>=Math.max(identityLower,30f),
					"Final identity acceptance must preserve the30percent floor and member lower: "+id);
			}else{
				require(identityLower==0 && blosumLower==0 && identityUpper==0 && blosumFinal==0,
					"Unavailable alignment gates must be explicit inert zeros: "+id);
			}
			float squared=0; double sum=0;
			for(float value:spectrum){
				if(!Float.isFinite(value) || value<0 || value>1){
					throw new IllegalArgumentException("Invalid spectrum frequency: "+id);
				}
				sum+=value; squared+=value*value;
			}
			// validWords counts K5 sketch words, not dimers: short canonical runs can have dimers but no K5 words.
			require(sum==0 || Math.abs(sum-1)<=0.000001,"Spectrum frequencies must sum to zero or one: "+id);
			spectrumNorm=(float)Math.sqrt(squared);
			assert(Float.isFinite(spectrumNorm) && (sum==0 || spectrumNorm>0)) : "Cached norm must support finite cosine: "+id;
		}
	}

	private static boolean upper(float value,float max){return Float.isFinite(value) && value>=0 && value<=max;}
	private static byte[] raw(ProteinSequence sequence){
		if(sequence==null){return null;}
		final byte[] result=new byte[sequence.enc.length];
		for(int i=0;i<result.length;i++){
			final byte code=sequence.enc[i];
			require(Blosum62.isStandard(code) || code==Blosum62.X_CODE,"Invalid encoded consensus residue");
			result[i]=code==Blosum62.X_CODE ? (byte)'X' : AminoAcid.numberToAcid[code];
		}
		return result;
	}
	private static void require(boolean condition,String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	public final String id;
	public final int rank;
	public final int scalarMask;
	public final double[] scalarBounds;
	public enum AlignmentCalibration{
		ENABLED("alignment_enabled"),
		UNAVAILABLE("alignment_unavailable");
		AlignmentCalibration(String label_){label=label_;}
		public final String label;
	}
	/** Explicit state is serialized only by v11/v12 live payloads. */
	public final AlignmentCalibration alignmentCalibration;
	final boolean explicitAlignmentCalibration;
	public final long members,residues,validWords;
	public final ProteinSequence consensus;
	public final byte[] consensusRaw;
	public final char[] ddl;
	/** 400 dimers or 20 exact projected row sums, as declared by the package header. */
	public final float[] spectrum;
	public final float spectrumNorm,cosineLower,identityLower,blosumLower,cosineAccept,identityUpper,blosumFinal;
}
