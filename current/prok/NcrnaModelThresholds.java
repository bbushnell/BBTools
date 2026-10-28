package prok;

import java.util.ArrayList;
import java.util.Arrays;
import fileIO.ByteFile;
import map.ObjectIntMap;
import map.ObjectMap;
import parse.LineParser1;

/** Immutable, name-bound identity and optional length cutoffs for RNA verification.
 * Every model in a named family is explicit; omitted families keep scalar defaults.
 * @author Raiden
 */
public final class NcrnaModelThresholds {
	public NcrnaModelThresholds(String[] names_,float[] pass_,float[] borderline_){
		this(names_,pass_,borderline_,null,null);
	}
	private NcrnaModelThresholds(String[] names_,float[] pass_,float[] borderline_,int[] minimum_,int[] maximum_){
		require(names_!=null && pass_!=null && borderline_!=null && names_.length>0 && names_.length==pass_.length && names_.length==borderline_.length,"Each model requires both identity cutoffs");
		names=names_.clone();pass=pass_.clone();borderline=borderline_.clone();final map.ObjectSet<String> seen=new map.ObjectSet<String>(String.class);
		require(minimum_==null && maximum_==null || minimum_!=null && maximum_!=null
			&& minimum_.length==names.length && maximum_.length==names.length,"Length cutoffs must bind the complete model roster");
		minimum=minimum_==null ? null : minimum_.clone();maximum=maximum_==null ? null : maximum_.clone();
		for(int i=0;i<names.length;i++){
			require(names[i]!=null && !names[i].isEmpty() && seen.add(names[i]),"Cutoffs must bind unique actual model names");
			require(Float.isFinite(pass[i]) && Float.isFinite(borderline[i]) && borderline[i]>=0 && pass[i]<=1 && borderline[i]<=pass[i],"Require0<=borderline<=pass<=1 for model "+names[i]);
			require(minimum==null || minimum[i]>=0 && maximum[i]>=Math.max(1,minimum[i]),"Model length bounds must be ordered and nonnegative: "+names[i]);
		}
	}
	/** Load all rows and validate all families before changing any live bundle. */
	public static void load(String path,ArrayList<NcrnaFamily> families){
		require(path!=null && families!=null,"Explicit cutoff table and loaded families required");
		final ObjectMap<String,Pending> pending=new ObjectMap<String,Pending>(String.class,Pending.class);final ArrayList<Pending> order=new ArrayList<Pending>();
		final ObjectMap<String,NcrnaFamily> known=new ObjectMap<String,NcrnaFamily>(String.class,NcrnaFamily.class);
		for(NcrnaFamily f:families){require(!known.contains(f.name),"Family names must uniquely bind model thresholds");known.put(f.name,f);}
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');int rows=0;
		try{RrnaResourceIO.header(in,"family\tmodel\tidpass\tidborderline",path);
			for(byte[] row;(row=in.nextLine())!=null;){p.set(row);require(p.terms()==4,"Cutoff table requires family/model/idpass/idborderline");
				final String name=p.parseString(0),model=p.parseString(1);final NcrnaFamily f=known.get(name);require(f!=null && f.reuseConsensusAlignment,"Model cutoffs require a loaded single-alignment family: "+name);
				Pending q=pending.get(name);if(q==null){q=new Pending(f);pending.put(name,q);order.add(q);}
				require(q.index.contains(model),"Unknown model in cutoff table: "+name+"/"+model);final int i=q.index.get(model);
				require(!q.seen[i],"Repeated model cutoff: "+name+"/"+model);q.seen[i]=true;
				q.pass[i]=(float)RrnaPositionalKmerTable.savedDouble(p,2);q.borderline[i]=(float)RrnaPositionalKmerTable.savedDouble(p,3);rows++;
			}
		}finally{require(!in.close(),"Model cutoff table read failed");}
		require(rows>0,"An empty override table cannot define a selected operating point");
		for(Pending q:order){for(boolean seen:q.seen){require(seen,"Every model in a named family requires an explicit row: "+q.family.name);}
			final NcrnaModelThresholds old=q.family.modelThresholds;
			if(old!=null){old.validate(q.family.modelNames);}
			q.result=new NcrnaModelThresholds(q.family.modelNames,q.pass,q.borderline,old==null ? null : old.minimum,old==null ? null : old.maximum);}
		for(Pending q:order){q.family.modelThresholds=q.result;for(int i=0;i<q.pass.length;i++){System.err.println("ncRNA model cutoff: family="+q.family.name+" model="+q.result.names[i]+" idPass="+q.pass[i]+" idBorderline="+q.borderline[i]);}}
	}
	/** General libraries bind every model. The historical mitochondrial coverage
	 * format binds added models only; original models retain unrestricted lengths.
	 * Bounds are already expanded by the resource producer, never scaled twice.
	 * An optional identity column overrides direct acceptance for that model.
	 * Its rescue floor is the smaller of the inherited floor and direct cutoff,
	 * so lowering direct acceptance cannot leave an unreachable passing interval. */
	public static void loadLengths(final String path, final NcrnaFamily family){
		require(path!=null && family!=null && family.modelNames!=null && family.library!=null
			&& family.modelNames.length==family.library.length, "Length profiles require a complete loaded model roster");
		final boolean coverage=MitoSsuModelPolicy.coverage(family.name,family.modelNames);
		require(!coverage || family.reuseConsensusAlignment, "Mitochondrial coverage retains its single-alignment verifier");
		final String[] names=family.modelNames;final int n=names.length,first=coverage ? NcrnaFamilyRegistry.get("16S").models.length : 0;
		final int[] minimum=new int[n],maximum=new int[n];Arrays.fill(maximum,Integer.MAX_VALUE);
		final float[] pass=new float[n],borderline=new float[n];
		final NcrnaModelThresholds old=family.modelThresholds;if(old!=null){old.validate(names);}
		for(int i=0; i<n; i++){pass[i]=old==null ? family.idPass : old.pass(i);borderline[i]=old==null ? family.idBorderline : old.borderline(i);}
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');int model=first;
		final boolean identities;
		try{
			final byte[] header=in.nextLine();
			require(header!=null, "A model profile requires its named column header: "+path);
			p.set(header);
			identities=p.terms()==5;
			require((p.terms()==4 || identities) && p.termEquals("model", 0)
				&& p.termEquals("model_length", 1) && p.termEquals("min_length", 2)
				&& p.termEquals("max_length", 3) && (!identities || p.termEquals("identity", 4)),
				"Expected model/model_length/min_length/max_length and optional identity: "+path);
			for(byte[] row;(row=in.nextLine())!=null;model++){
				p.set(row);require(p.terms()==(identities ? 5 : 4) && model<n,"Model profile must contain exactly one complete row per selected model: "+path);
				require(p.termEquals(names[model],0) && p.parseInt(1)==family.library[model].length,
					"Length profile name/order/sequence length differs from caller model "+model+": "+path);
				minimum[model]=p.parseInt(2);maximum[model]=p.parseInt(3);
				require(minimum[model]>0 && maximum[model]>=minimum[model],"Positive ordered gene-length limits required for "+names[model]);
				if(identities){
					pass[model]=(float)RrnaPositionalKmerTable.savedDouble(p, 4);
					require(Float.isFinite(pass[model]) && pass[model]>=0 && pass[model]<=1,
						"Model acceptance identity must be finite in [0,1]: "+names[model]);
					borderline[model]=Math.min(borderline[model], pass[model]);
				}
			}
		}finally{require(!in.close(),"Model length profile read failed: "+path);}
		require(model==n,"Missing model length rows: "+path);
		family.modelThresholds=new NcrnaModelThresholds(names,pass,borderline,minimum,maximum);
		if(identities){
			float low=1, high=0;
			for(int i=first; i<n; i++){low=Math.min(low, pass[i]);high=Math.max(high, pass[i]);}
			System.err.println("ncRNA model identities: family="+family.name+" models="+(n-first)+" min="+low+" max="+high+" profile="+path);
		}
	}
	void validate(String[] actual){require(Arrays.equals(names,actual),"Model cutoff order must exactly match the caller consensus library");}
	static final class Pending{
		Pending(NcrnaFamily f){family=f;require(f.modelNames!=null && f.library!=null && f.modelNames.length==f.library.length,"Loaded family must bind every consensus to a model name");final int n=f.modelNames.length;pass=new float[n];borderline=new float[n];seen=new boolean[n];
			for(int i=0;i<n;i++){require(!index.contains(f.modelNames[i]),"Duplicate caller model name");index.put(f.modelNames[i],i);}}
		final NcrnaFamily family;final ObjectIntMap<String> index=new ObjectIntMap<String>(String.class);final float[] pass,borderline;final boolean[] seen;NcrnaModelThresholds result;
	}
	static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	private final String[] names;
	private final float[] pass,borderline;
	private final int[] minimum,maximum;
	public float pass(int model){return pass[model];}
	public float borderline(int model){return borderline[model];}
	public int minimum(int model){return minimum==null ? 0 : minimum[model];}
	public int maximum(int model){return maximum==null ? Integer.MAX_VALUE : maximum[model];}
	public boolean hasLengths(){return minimum!=null;}
	/** A four-column length profile retains inherited identity thresholds exactly. */
	public boolean sameIdentities(float scalarPass, float scalarBorderline){
		for(int i=0; i<pass.length; i++){if(pass[i]!=scalarPass || borderline[i]!=scalarBorderline){return false;}}return true;
	}
}
