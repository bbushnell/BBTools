package prok;

import java.util.ArrayList;
import java.util.Arrays;
import fileIO.ByteFile;
import map.ObjectIntMap;
import map.ObjectMap;
import parse.LineParser1;

/** Immutable, name-bound identity cutoffs for the single-alignment RNA path.
 * Every model in a named family is explicit; omitted families keep scalar defaults.
 * @author Raiden
 */
public final class NcrnaModelThresholds {
	NcrnaModelThresholds(String[] names_,float[] pass_,float[] borderline_){
		require(names_!=null && pass_!=null && borderline_!=null && names_.length>0 && names_.length==pass_.length && names_.length==borderline_.length,"Each model requires both identity cutoffs");
		names=names_.clone();pass=pass_.clone();borderline=borderline_.clone();final map.ObjectSet<String> seen=new map.ObjectSet<String>(String.class);
		for(int i=0;i<names.length;i++){
			require(names[i]!=null && !names[i].isEmpty() && seen.add(names[i]),"Cutoffs must bind unique actual model names");
			require(Float.isFinite(pass[i]) && Float.isFinite(borderline[i]) && borderline[i]>=0 && pass[i]<=1 && borderline[i]<=pass[i],"Require0<=borderline<=pass<=1 for model "+names[i]);
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
		for(Pending q:order){for(boolean seen:q.seen){require(seen,"Every model in a named family requires an explicit row: "+q.family.name);}q.result=new NcrnaModelThresholds(q.family.modelNames,q.pass,q.borderline);}
		for(Pending q:order){q.family.modelThresholds=q.result;for(int i=0;i<q.pass.length;i++){System.err.println("ncRNA model cutoff: family="+q.family.name+" model="+q.result.names[i]+" idPass="+q.pass[i]+" idBorderline="+q.borderline[i]);}}
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
	float pass(int model){return pass[model];}
	float borderline(int model){return borderline[model];}
}
