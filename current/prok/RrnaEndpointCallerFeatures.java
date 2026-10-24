package prok;

import java.util.Arrays;
import idaligner.AlignmentStats;
import ml.CellNet;

/** Per-scavenger scratch for the measured rRNA endpoint feature contract.
 * Coordinates are inclusive in the current strand-oriented contig. Extraction
 * observes verifier-accepted candidates before snapshot/DP filtering. Optional
 * explicitly supplied networks refine ends; observation alone never moves them.
 * @author Raiden
 */
public final class RrnaEndpointCallerFeatures {
	public static final class Resources {
		/** Preserve the public euk5S API while new families supply their identity explicitly. */
		public Resources(String[] names_,byte[][] refs_,RrnaPositionalKmerTable.Table[] five_,RrnaPositionalKmerTable.Table[] three_){
			this(legacyFamily(names_),names_,refs_,five_,three_,null,null);
		}
		public Resources(String[] names_,byte[][] refs_,RrnaPositionalKmerTable.Table[] five_,RrnaPositionalKmerTable.Table[] three_,CellNet[] fiveNets_,CellNet[] threeNets_){
			this(legacyFamily(names_),names_,refs_,five_,three_,fiveNets_,threeNets_);
		}
		static String legacyFamily(String[] names){
			require(names!=null && names.length>0,"Legacy endpoint resources bind a nonempty euk5S library");
			for(String name:names){require(name!=null && name.startsWith("euk5S_"),"Legacy endpoint resources require euk5S consensus names");}
			return "euk5S";
		}
		public Resources(String family_,String[] names_,byte[][] refs_,RrnaPositionalKmerTable.Table[] five_,RrnaPositionalKmerTable.Table[] three_){
			this(family_,names_,refs_,five_,three_,null,null);
		}
		public Resources(String family_,String[] names_,byte[][] refs_,RrnaPositionalKmerTable.Table[] five_,RrnaPositionalKmerTable.Table[] three_,CellNet[] fiveNets_,CellNet[] threeNets_){
			this(family_,names_,refs_,five_,three_,fiveNets_,threeNets_,null);
		}
		public Resources(String family_,String[] names_,byte[][] refs_,RrnaPositionalKmerTable.Table[] five_,RrnaPositionalKmerTable.Table[] three_,CellNet[] fiveNets_,CellNet[] threeNets_,boolean[] raw_){
			this(family_,names_,refs_,five_,three_,fiveNets_,threeNets_,raw_,false);
		}
		public Resources(String family_,String[] names_,byte[][] refs_,RrnaPositionalKmerTable.Table[] five_,RrnaPositionalKmerTable.Table[] three_,CellNet[] fiveNets_,CellNet[] threeNets_,boolean[] raw_,boolean acceptedIdentity_){
			require(family_!=null && family_.matches("[A-Za-z0-9_]+"),"Endpoint resources require an explicit registered family key");family=family_;
			require(!acceptedIdentity_ || supportsAcceptedIdentity(family),"Only explicitly trained16S/23S/prok5S resources may select accepted-verifier identity; other families keep their candidate-span definition");acceptedVerifierIdentity=acceptedIdentity_;
			require(names_!=null && refs_!=null && five_!=null && three_!=null && names_.length>0 && names_.length==refs_.length && names_.length==five_.length && names_.length==three_.length,"Every caller consensus must bind one table per end");
			names=names_.clone();refs=refs_.clone();five=five_.clone();three=three_.clone();final map.ObjectSet<String> seen=new map.ObjectSet<String>(String.class);
			require(raw_==null || raw_.length==names.length,"Explicit raw-endpoint mask must bind every detection model");raw=raw_==null?new boolean[names.length]:raw_.clone();
			require((fiveNets_==null && threeNets_==null) || (fiveNets_!=null && threeNets_!=null && fiveNets_.length==names.length && threeNets_.length==names.length),"Inference requires both networks for every model");
			fiveNets=fiveNets_==null?null:fiveNets_.clone();threeNets=threeNets_==null?null:threeNets_.clone();
			int first=-1;for(int i=0;i<names.length;i++){if(!raw[i]){first=i;break;}}
			require(first>=0 && five[first]!=null,"At least one refined model must declare the family's endpoint geometry");
			radius=five[first].radius;sites=five[first].sites;inputs=five[first].inputs;
			for(int i=0;i<names.length;i++){
				require(names[i]!=null && names[i].matches("[A-Za-z0-9_]+") && seen.add(names[i]) && refs[i]!=null && refs[i].length>0,"Feature resources require unique safe caller model names and actual consensus bases");
				if(raw[i]){require(five[i]==null && three[i]==null && (fiveNets==null || (fiveNets[i]==null && threeNets[i]==null)),"Explicit raw models must not silently load ignored tables or networks");continue;}
				require(five[i]!=null && three[i]!=null && five[i].model.equals(names[i]) && three[i].model.equals(names[i]) && five[i].end.equals("5prime") && three[i].end.equals("3prime"),"Each table must bind the trained model and its own endpoint");
				// Geometry is shared within a family; k and anchor remain per end/model.
				require(five[i].radius==radius && three[i].radius==radius,"All model/end tables in one family must use its declared radius");
				require(five[i].finished && three[i].finished,"Shared score tables must be finalized before caller workers start");
				if(fiveNets!=null){validateNet(fiveNets[i],inputs,sites);validateNet(threeNets[i],inputs,sites);}
			}
		}
		final String family;public final int radius,sites,inputs;public final boolean acceptedVerifierIdentity;
		public final boolean[] raw;
		public final String[] names;public final byte[][] refs;public final RrnaPositionalKmerTable.Table[] five,three;
		public final CellNet[] fiveNets,threeNets;
		static void validateNet(CellNet net,int inputs,int sites){require(net!=null && net.numInputs()==inputs && net.numOutputs()==sites,"Endpoint network dimensions must match the family's table radius and three extras");}
		static void validateNet(CellNet net){validateNet(net,28,25);}
		static boolean supportsAcceptedIdentity(String family){return "16S".equals(family) || "23S".equals(family) || "prok5S".equals(family);}
	}
	/** Called synchronously with reusable worker-local arrays. Copy before retaining;
	 * implementations shared by multiple workers must synchronize their own output. */
	public interface Sink {
		void capture(String contig,int strand,int model,int rawStart,int rawStop,
			float[] five,boolean fiveUsable,float[] three,boolean threeUsable);
	}
	public RrnaEndpointCallerFeatures(Resources resources_,Sink sink_){
		require(resources_!=null && (sink_!=null || resources_.fiveNets!=null),"Endpoint processing requires explicit observation or inference");resources=resources_;sink=sink_;
		fiveNets=copyNets(resources.fiveNets);threeNets=copyNets(resources.threeNets);
		scratch=new float[resources.sites];five=new float[resources.inputs];three=new float[resources.inputs];
	}
	public void capture(Orf candidate,byte[] bases,int model,float contigGC){
		capture(candidate,bases,model,contigGC,1,Integer.MAX_VALUE);
	}
	public void capture(Orf candidate,byte[] bases,int model,float contigGC,int minLen,int maxLen){
		capture(candidate,bases,model,contigGC,minLen,maxLen,Float.NaN);
	}
	public void capture(Orf candidate,byte[] bases,int model,float contigGC,int minLen,int maxLen,float acceptedIdentity){
		capture(candidate,bases,model,contigGC,minLen,maxLen,acceptedIdentity,0);
	}
	/** Fixed endpoints are in the current strand-oriented frame: bit1 is5prime,
	 * bit2 is3prime. Disabled ends are neither scored nor predicted. The caller
	 * supplies its already-established partial length floor; training reference
	 * length and accepted-verifier identity retain their original definitions. */
	public void capture(Orf candidate,byte[] bases,int model,float contigGC,int minLen,int maxLen,float acceptedIdentity,int fixedMask){
		require((fixedMask&~3)==0,"Endpoint fixed mask permits only oriented5prime/3prime bits");
		require(candidate!=null && resources.family.equals(candidate.ncrnaFamily) && model>=0 && model<resources.refs.length && resources.names[model].equals(candidate.trnaModel),"Raw candidate must belong to the configured family and winning model");
		if(fixedMask==3){return;}// No endpoint features, identity realignment, inference or observation.
		final byte[] reference=resources.refs[model];final float identity;
		if(resources.acceptedVerifierIdentity){require(Float.isFinite(acceptedIdentity) && acceptedIdentity>=0 && acceptedIdentity<=1,"Accepted-identity features require the actual winning-model verifier identity; no cropped re-alignment fallback");identity=acceptedIdentity;}
		else{identity=RrnaEndpointVector.candidateIdentity(bases,candidate.start,candidate.stop,reference,stats);}
		// Legacy bundles retain their counted candidate-span alignment. Accepted-
		// identity bundles reuse the supplied scalar; neither path invents endpoint
		// predictions for explicit raw models.
		if(resources.raw[model]){return;}
		if((fixedMask&1)!=0){Arrays.fill(five,Float.NaN);}
		if((fixedMask&2)!=0){Arrays.fill(three,Float.NaN);}
		final boolean fiveUsable=(fixedMask&1)==0 && NcrnaBoundaryScorer.rrnaEndpointFeatures(resources.five[model],bases,candidate.start,candidate.stop,reference.length,contigGC,identity,scratch,five);
		final boolean threeUsable=(fixedMask&2)==0 && NcrnaBoundaryScorer.rrnaEndpointFeatures(resources.three[model],bases,candidate.start,candidate.stop,reference.length,contigGC,identity,scratch,three);
		if(sink!=null){sink.capture(candidate.scafName,candidate.strand,model,candidate.start,candidate.stop,five,fiveUsable,three,threeUsable);}
		if(fiveNets!=null){
			final int start=fiveUsable?RrnaEndpointVector.endpointFromClass(candidate.start,predict(fiveNets[model],five),resources.radius):candidate.start;
			final int stop=threeUsable?RrnaEndpointVector.endpointFromClass(candidate.stop,predict(threeNets[model],three),resources.radius):candidate.stop;
			apply(candidate,start,stop,bases.length,minLen,maxLen);
		}
	}
	/** Independent argmax at each raw end; a jointly invalid span keeps both raw ends.
	 * The pre-refinement detection score remains unchanged, as in the older NN path. */
	public static boolean apply(Orf candidate,int start,int stop,int length,int minLen,int maxLen){
		assert(minLen>0 && maxLen>=minLen):shared.KillSwitch.assertDie("Endpoint fallback uses the caller's valid length interval");
		final long span=stop-start+1L;
		if(start<0 || stop>=length || span<minLen || span>maxLen){return false;}
		candidate.start=start;candidate.stop=stop;return true;
	}
	public static int predict(CellNet net,float[] features){net.applyInput(features);net.feedForward();return RrnaResourceIO.argmax(net.getOutput());}
	static CellNet[] copyNets(CellNet[] templates){if(templates==null){return null;}final CellNet[] copies=new CellNet[templates.length];for(int i=0;i<copies.length;i++){copies[i]=templates[i]==null?null:templates[i].copy(false);}return copies;}
	static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	final Resources resources;final Sink sink;final AlignmentStats stats=new AlignmentStats(true);
	public final CellNet[] fiveNets,threeNets;
	public final float[] scratch,five,three;
}
