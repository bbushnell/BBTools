package prok;

/** Shared member voting for plain DDL geometries with optional pooled consensus
 * bonuses. One word definition enumerates each real member once; layouts with the
 * same member bins share their transient vote map across bonus choices. Finished
 * vote selection yields native keys for wrapping in DDL reference objects.
 * Native insertion owns all hashing/binning. No support/origin arrays reach a scorer.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlLayoutModels{

	/*--------------------------------------------------------------*/
	/*----------------          Construction        ----------------*/
	/*--------------------------------------------------------------*/

	public static Result buildMembers(final NcrnaDdlMembers members, final NcrnaDdlLayoutKeys[] layouts){
		if(layouts==null){throw new IllegalArgumentException("Member voting requires an explicit layout roster");}
		return buildMembers(members, layouts, null, new int[layouts.length], new boolean[layouts.length]);
	}

	public static Result buildMembers(final NcrnaDdlMembers members, final NcrnaDdlLayoutKeys[] layouts,
			final byte[][] consensuses, final int[] bonuses, final boolean[] consensusShort){
		if(members==null || members.models()<1 || layouts==null || layouts.length<1
				|| bonuses==null || consensusShort==null || bonuses.length!=layouts.length || consensusShort.length!=layouts.length){
			throw new IllegalArgumentException("Member voting requires actual models and aligned geometry/bonus declarations");
		}
		final int[] owners=new int[layouts.length];
		final NcrnaDdlLayoutKeys[] bonusKeys=new NcrnaDdlLayoutKeys[layouts.length];
		boolean needConsensus=false;
		for(int geometry=0; geometry<layouts.length; geometry++){
			final NcrnaDdlLayoutKeys layout=layouts[geometry];
			if(layout==null || bonuses[geometry]<0){
				throw new IllegalArgumentException("Every member geometry needs an explicit codec and nonnegative consensus bonus");
			}
			owners[geometry]=geometry;
			for(int before=0; before<geometry; before++){
				if(layout.sameGeometry(layouts[before])){
					owners[geometry]=owners[before];
					break;
				}
			}
			if(bonuses[geometry]>0){
				needConsensus=true;
				if(consensuses==null || consensuses.length!=members.models()){
					throw new IllegalArgumentException("A pooled bonus must bind every actual model consensus");
				}
				bonusKeys[geometry]=new NcrnaDdlLayoutKeys(layout.bins(), layout.shortBins(), consensusShort[geometry],
					layout.format(), layout.shortK(), layout.longK(), layout.maskedLong());
			}
		}
		final NcrnaDdlLayoutKeys first=layouts[0];
		final NcrnaDdlLayoutKeys.Raw raw=new NcrnaDdlLayoutKeys.Raw(first.shortK(), first.longK(), first.maskedLong());
		for(final NcrnaDdlLayoutKeys layout:layouts){
			if(layout.shortK()!=raw.shortK || layout.longK()!=raw.longK || layout.maskedLong()!=raw.maskedLong){
				throw new IllegalArgumentException("One shared pass needs one exact short/long/mask definition");
			}
		}
		final NcrnaDdlLayoutVotes.Result[][] profiles=new NcrnaDdlLayoutVotes.Result[layouts.length][members.models()];
		long voters=0;
		int shortExponent=-1, longExponent=-1;
		for(int model=0; model<members.models(); model++){
			final byte[][] sequences=members.members(model);
			if(sequences==null || sequences.length<1){throw new IllegalArgumentException("Every model needs actual retained voters: model="+model);}
			final NcrnaDdlLayoutVotes[] votes=new NcrnaDdlLayoutVotes[layouts.length];
			for(int geometry=0; geometry<layouts.length; geometry++){
				if(owners[geometry]==geometry){votes[geometry]=new NcrnaDdlLayoutVotes(layouts[geometry].bins());}
			}
			for(final byte[] sequence:sequences){
				raw.build(sequence);
				shortExponent=Math.max(shortExponent, raw.maxExponent(true));
				longExponent=Math.max(longExponent, raw.maxExponent(false));
				voters++;
				for(int geometry=0; geometry<layouts.length; geometry++){
					if(owners[geometry]!=geometry){continue;}
					layouts[geometry].build(raw);
					votes[geometry].add(layouts[geometry].maximumKeys());
				}
			}
			if(needConsensus){raw.build(consensuses[model]);}
			for(int geometry=0; geometry<layouts.length; geometry++){
				final NcrnaDdlLayoutKeys consensus=bonusKeys[geometry];
				if(consensus!=null){consensus.build(raw);}
				profiles[geometry][model]=votes[owners[geometry]].finish(consensus, bonuses[geometry]);
				assert(profiles[geometry][model].members==sequences.length)
					: "Every bonus choice must preserve the same real voters after shared hashing and voting";
			}
		}
		assert(voters==members.keptMembers()) : "Each kept member must be hashed once per declared word definition";
		return new Result(profiles, voters, shortExponent, longExponent);
	}

	/*--------------------------------------------------------------*/
	/*----------------             Result           ----------------*/
	/*--------------------------------------------------------------*/

	public static final class Result{
		private Result(final NcrnaDdlLayoutVotes.Result[][] profiles_, final long voters_, final int short_, final int long_){
			profiles=profiles_;
			voters=voters_;
			shortExponent=short_;
			longExponent=long_;
		}
		public final NcrnaDdlLayoutVotes.Result[][] profiles;
		public final long voters;
		public final int shortExponent, longExponent;
	}

	private NcrnaDdlLayoutModels(){}
}
