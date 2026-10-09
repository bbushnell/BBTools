package prok;

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
			require(family_!=null && family_.matches("[A-Za-z0-9_]+"),"Endpoint resources require an explicit registered family key");family=family_;
			require(names_!=null && refs_!=null && five_!=null && three_!=null && names_.length>0 && names_.length==refs_.length && names_.length==five_.length && names_.length==three_.length,"Every caller consensus must bind one table per end");
			names=names_.clone();refs=refs_.clone();five=five_.clone();three=three_.clone();final map.ObjectSet<String> seen=new map.ObjectSet<String>(String.class);
			require((fiveNets_==null && threeNets_==null) || (fiveNets_!=null && threeNets_!=null && fiveNets_.length==names.length && threeNets_.length==names.length),"Inference requires both networks for every model");
			fiveNets=fiveNets_==null?null:fiveNets_.clone();threeNets=threeNets_==null?null:threeNets_.clone();
			require(five[0]!=null,"The first model must declare the family's endpoint geometry");
			radius=five[0].radius;sites=five[0].sites;inputs=five[0].inputs;
			for(int i=0;i<names.length;i++){
				require(names[i]!=null && names[i].matches("[A-Za-z0-9_]+") && seen.add(names[i]) && refs[i]!=null && refs[i].length>0,"Feature resources require unique safe caller model names and actual consensus bases");
				require(five[i]!=null && three[i]!=null && five[i].model.equals(names[i]) && three[i].model.equals(names[i]) && five[i].end.equals("5prime") && three[i].end.equals("3prime"),"Each table must bind the trained model and its own endpoint");
				// Geometry is shared within a family; k and anchor remain per end/model.
				require(five[i].radius==radius && three[i].radius==radius,"All model/end tables in one family must use its declared radius");
				require(five[i].finished && three[i].finished,"Shared score tables must be finalized before caller workers start");
				if(fiveNets!=null){validateNet(fiveNets[i],inputs,sites);validateNet(threeNets[i],inputs,sites);}
			}
		}
		final String family;final int radius,sites,inputs;
		final String[] names;final byte[][] refs;final RrnaPositionalKmerTable.Table[] five,three;
		final CellNet[] fiveNets,threeNets;
		static void validateNet(CellNet net,int inputs,int sites){require(net!=null && net.numInputs()==inputs && net.numOutputs()==sites,"Endpoint network dimensions must match the family's table radius and three extras");}
		static void validateNet(CellNet net){validateNet(net,28,25);}
	}
	/** Called synchronously with reusable worker-local arrays. Copy before retaining;
	 * implementations shared by multiple workers must synchronize their own output. */
	public interface Sink {
		void capture(String contig,int strand,int model,int rawStart,int rawStop,
			float[] five,boolean fiveUsable,float[] three,boolean threeUsable);
	}
	RrnaEndpointCallerFeatures(Resources resources_,Sink sink_){
		require(resources_!=null && (sink_!=null || resources_.fiveNets!=null),"Endpoint processing requires explicit observation or inference");resources=resources_;sink=sink_;
		fiveNets=copyNets(resources.fiveNets);threeNets=copyNets(resources.threeNets);
		scratch=new float[resources.sites];five=new float[resources.inputs];three=new float[resources.inputs];
	}
	void capture(Orf candidate,byte[] bases,int model,float contigGC){
		capture(candidate,bases,model,contigGC,1,Integer.MAX_VALUE);
	}
	void capture(Orf candidate,byte[] bases,int model,float contigGC,int minLen,int maxLen){
		require(candidate!=null && resources.family.equals(candidate.ncrnaFamily) && model>=0 && model<resources.refs.length && resources.names[model].equals(candidate.trnaModel),"Raw candidate must belong to the configured family and winning model");
		final byte[] reference=resources.refs[model];final float identity=RrnaEndpointVector.candidateIdentity(bases,candidate.start,candidate.stop,reference,stats);
		final boolean fiveUsable=NcrnaBoundaryScorer.rrnaEndpointFeatures(resources.five[model],bases,candidate.start,candidate.stop,reference.length,contigGC,identity,scratch,five);
		final boolean threeUsable=NcrnaBoundaryScorer.rrnaEndpointFeatures(resources.three[model],bases,candidate.start,candidate.stop,reference.length,contigGC,identity,scratch,three);
		if(sink!=null){sink.capture(candidate.scafName,candidate.strand,model,candidate.start,candidate.stop,five,fiveUsable,three,threeUsable);}
		if(fiveNets!=null){
			final int start=fiveUsable?RrnaEndpointVector.endpointFromClass(candidate.start,predict(fiveNets[model],five),resources.radius):candidate.start;
			final int stop=threeUsable?RrnaEndpointVector.endpointFromClass(candidate.stop,predict(threeNets[model],three),resources.radius):candidate.stop;
			apply(candidate,start,stop,bases.length,minLen,maxLen);
		}
	}
	/** Independent argmax at each raw end; a jointly invalid span keeps both raw ends.
	 * The pre-refinement detection score remains unchanged, as in the older NN path. */
	static boolean apply(Orf candidate,int start,int stop,int length,int minLen,int maxLen){
		assert(minLen>0 && maxLen>=minLen):shared.KillSwitch.assertDie("Endpoint fallback uses the caller's valid length interval");
		final long span=stop-start+1L;
		if(start<0 || stop>=length || span<minLen || span>maxLen){return false;}
		candidate.start=start;candidate.stop=stop;return true;
	}
	static int predict(CellNet net,float[] features){net.applyInput(features);net.feedForward();return RrnaResourceIO.argmax(net.getOutput());}
	static CellNet[] copyNets(CellNet[] templates){if(templates==null){return null;}final CellNet[] copies=new CellNet[templates.length];for(int i=0;i<copies.length;i++){copies[i]=templates[i].copy(false);}return copies;}
	static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	final Resources resources;final Sink sink;final AlignmentStats stats=new AlignmentStats(true);
	final CellNet[] fiveNets,threeNets;
	final float[] scratch,five,three;
}
