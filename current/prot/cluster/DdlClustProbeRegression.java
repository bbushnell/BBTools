package prot.cluster;

import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntHashMap;

/** Fits local structural-objective scales from isolated matched-control probes.
 * This joins raw proposal features with compact reassignment responses only;
 * no sequence data, family codecs, coefficient search, or candidate selection
 * occurs here.
 * @author Yelan
 */
public final class DdlClustProbeRegression{

	public static void main(final String[] args) throws Exception{
		final DdlClustProbeRegression tool=new DdlClustProbeRegression(args);
		try{tool.process();}finally{Shared.closeStream(tool.outstream);}
	}

	public DdlClustProbeRegression(String[] args){
		final PreParser pre=new PreParser(args,getClass(),false);
		args=pre.args; outstream=pre.outstream;
		final Parser parser=new Parser(); parser.overwrite=false;
		for(String arg:args){
			final String[] pair=Parse.splitOnFirst(arg,'=');
			final String key=pair[0].toLowerCase(Locale.ROOT);
			final String value=pair.length>1 && !pair[1].equalsIgnoreCase("null") ? pair[1] : null;
			if(key.equals("results")){results=value;}
			else if(key.equals("report")){report=value;}
			else if(key.equals("localresults")){localResults=value;}
			else if(key.equals("featuretable")){featureTable=value;}
			else if(key.equals("fitslot")){fitSlot=Parse.parseBoolean(value);}
			else if(key.equals("fitsize")){fitSize=Parse.parseBoolean(value);}
			else if(key.equals("outresiduals")){outResiduals=value;}
			else if(!parser.parse(arg,key,value)){throw new IllegalArgumentException("Unknown probe-regression parameter: "+arg);}
		}
		out=parser.out1; overwrite=parser.overwrite;
		require(!parser.append,"Probe regressions cannot append");
		require(out!=null,"Require out");
		require(featureTable!=null ? localResults==null && results==null && report==null :
			(localResults!=null && results==null && report==null) ||
			(localResults==null && results!=null && report!=null),
			"Require exactly one of featuretable, localresults, or legacy results plus report");
		require(featureTable==null || (!fitSize && !fitSlot),"featuretable requires fitsize=f and fitslot=f");
	}

	public void process() throws IOException{
		final ArrayList<Path> paths=new ArrayList<Path>();
		if(featureTable!=null){register(paths,featureTable,true);}
		else if(localResults==null){
			register(paths,results,true); register(paths,report,true);
		}else{register(paths,localResults,true);}
		register(paths,out,false);
		if(outResiduals!=null){register(paths,outResiduals,false);}
		require(overwrite || !Files.exists(Paths.get(out),LinkOption.NOFOLLOW_LINKS),"Regression output exists and overwrite=false: "+out);
		require(Tools.testOutputFiles(overwrite,false,false,out),"Probe-regression output is unavailable");
		if(outResiduals!=null){
			require(overwrite || !Files.exists(Paths.get(outResiduals),LinkOption.NOFOLLOW_LINKS),
				"Residual output exists and overwrite=false: "+outResiduals);
			require(Tools.testOutputFiles(overwrite,false,false,outResiduals),"Residual output is unavailable");
		}
		final ArrayList<Probe> probes=featureTable!=null ? readFeatureTable() :
			(localResults==null ? readResults() : readLocalResults());
		if(featureTable==null && localResults==null){readResponses(probes);}
		final ArrayList<Probe> observations=new ArrayList<Probe>();
		for(Probe p:probes){if(p.fitEligible){observations.add(p);}}
		final Fit fit=fit(observations,fitSize,fitSlot);
		write(fit,observations.size(),probes.size());
		if(outResiduals!=null){writeResiduals(probes);}
		outstream.println("DDLCLUST_PROBE_REGRESSION_PASS observations="+observations.size()+" predictors="+fit.coefficients.length);
	}

	private ArrayList<Probe> readResults(){
		final ArrayList<Probe> rows=new ArrayList<Probe>();
		final IntHashMap ids=new IntHashMap();
		final ByteFile reader=ByteFile.makeByteFile(results,false);
		try{
			header(reader.nextLine(),StructuralProbeBatch.RESULT_HEADER,"producer");
			final LineParser1 row=new LineParser1('\t');
			for(byte[] line=reader.nextLine(); line!=null; line=reader.nextLine()){
				row.set(line); require(row.terms()==31,"Producer rows require exactly 31 columns");
				final int id=identifier(row,0,false);
				require(ids.get(id)<0,"Duplicate candidate_index in producer results"); ids.put(id,rows.size());
				final boolean split=row.termEquals("fission",2);
				require(split || row.termEquals("fusion",2),"Unknown operation in producer results");
				final long feasible=nonnegative(row,3);
				require(feasible<=1,"Feasible must be zero or one in producer results");
				final int a=identifier(row,5,false), b=identifier(row,6,true);
				final long initial=nonnegative(row,30);
				require(initial>0,"Initial affected member count must be positive");
				final Probe p=new Probe(id,a,b,split,feasible==1,initial);
				if(p.feasible){
					p.cohesion=finite(row,20); p.separation=finite(row,21);
					p.size=finite(row,22); p.slot=finite(row,23);
				}else{
					for(int i=19; i<=23; i++){require(row.termEquals("NA",i),"Rejected producer row retained fitted predictors");}
				}
				rows.add(p);
			}
		}finally{close(reader,results);}
		require(!rows.isEmpty(),"Empty producer table");
		return rows;
	}

	private void readResponses(final ArrayList<Probe> probes){
		final IntHashMap ids=new IntHashMap(probes.size());
		for(int i=0; i<probes.size(); i++){ids.put(probes.get(i).id,i);}
		final ByteFile reader=ByteFile.makeByteFile(report,false);
		try{
			header(reader.nextLine(),DdlClustProbeReport.REPORT_HEADER,"probe report");
			final LineParser1 row=new LineParser1('\t');
			long inputRecords=-1;
			int seen=0;
			for(byte[] line=reader.nextLine(); line!=null; line=reader.nextLine()){
				row.set(line); require(row.terms()==13,"Probe-report rows require exactly 13 columns");
				final int id=identifier(row,0,false),offset=ids.get(id);
				require(offset>=0,"Probe report contains an unknown candidate_index");
				final Probe p=probes.get(offset);
				require(!p.seen,"Duplicate candidate_index in probe report"); p.seen=true; seen++;
				require(row.termEquals(p.split ? "fission" : "fusion",1),"Probe-report operation disagrees with producer");
				final long feasible=nonnegative(row,2);
				require(feasible==(p.feasible ? 1 : 0),"Probe-report feasibility disagrees with producer");
				require(identifier(row,4,false)==p.a && identifier(row,5,true)==p.b,
					"Probe-report endpoints disagree with producer");
				final long initial=nonnegative(row,6), control=nonnegative(row,7);
				final long input=nonnegative(row,12);
				require(initial==p.initial,"Probe-report initial count disagrees with producer");
				require(initial<=input && control<=input,"Probe-report counts exceed the input denominator");
				if(inputRecords<0){inputRecords=input;}
				else{require(input==inputRecords,"Probe-report input denominators disagree");}
				final long controlDelta=decimal(row,10);
				require(controlDelta==Math.subtractExact(control,initial),
					"Probe-report control-minus-initial arithmetic disagrees");
				if(p.feasible){
					final long modified=nonnegative(row,8), post=decimal(row,9), response=decimal(row,11);
					require(modified<=input,"Probe-report counts exceed the input denominator");
					require(post==Math.subtractExact(modified,initial) &&
						response==Math.subtractExact(modified,control),
						"Probe-report response arithmetic disagrees");
					p.response=(double)response;
				}else{
					require(row.termEquals("NA",8) && row.termEquals("NA",9) &&
						row.termEquals("NA",11),"Rejected probe has a fitted response");
				}
			}
			require(seen==probes.size(),"Probe report omitted producer rows");
		}finally{close(reader,report);}
	}

	private ArrayList<Probe> readLocalResults(){
		final ArrayList<Probe> rows=new ArrayList<Probe>();
		final IntHashMap ids=new IntHashMap();
		final ByteFile reader=ByteFile.makeByteFile(localResults,false);
		try{
			header(reader.nextLine(),LocalProbeBatch.HEADER,"local response");
			final LineParser1 row=new LineParser1('\t');
			for(byte[] line=reader.nextLine(); line!=null; line=reader.nextLine()){
				row.set(line); require(row.terms()==29,"Local response rows require exactly 29 columns");
				require(row.termEquals("local_replacements_v1",0),"Unsupported local response kind");
				final int id=identifier(row,1,false);
				final long seed=decimal(row,2);
				require(ids.get(id)<0,"Duplicate candidate_index in local response"); ids.put(id,rows.size());
				final boolean split=row.termEquals("fission",3);
				require(split || row.termEquals("fusion",3),"Unknown operation in local response");
				final int a=identifier(row,4,false), b=identifier(row,5,true);
				require(split ? b==-1 : b>=0 && b!=a,"Local response endpoints disagree with operation");
				final long feasible=nonnegative(row,6);
				require(feasible<=1,"Feasible must be zero or one in local response");
				final long affected=nonnegative(row,8), poolSize=nonnegative(row,9);
				require(affected>0,"Local response affected count must be positive");
				require(affected<=Integer.MAX_VALUE && poolSize<=Integer.MAX_VALUE &&
					poolSize+affected<=Integer.MAX_VALUE,
					"Local response pool plus affected count exceeds the native int input range");
				final long slotCost=nonnegative(row,10);
				require(slotCost<=Integer.MAX_VALUE,"Local response slot cost exceeds the native int range");
				final double temperature=finite(row,11);
				require(temperature>=0,"Local response temperature cannot be negative");
				nonnegative(row,28);
				final Probe p=new Probe(id,a,b,split,feasible==1,affected); p.seen=true; p.sampleSeed=Long.toString(seed);
				if(p.feasible){
					final double balance=finite(row,12);
					require(balance>=0 && balance<=1,"Local response balance is outside [0,1]");
					p.cohesion=finite(row,13); p.separation=finite(row,14);
					p.size=finite(row,15); p.slot=finite(row,16);
					require(p.slot==(split ? -slotCost : slotCost),
						"Local response slot delta disagrees with operation and slot cost");
					final long netGain=decimal(row,17);
					final long poolVisited=nonnegative(row,18), poolAccepted=nonnegative(row,19);
					final long affectedAccepted=nonnegative(row,20), accepted=nonnegative(row,21);
					final long candidates=nonnegative(row,22), alignments=nonnegative(row,23);
					final long lowHits=nonnegative(row,24), hydro=nonnegative(row,25);
					nonnegative(row,26);
					final long affectedLost=nonnegative(row,27);
					require(poolVisited<=poolSize && poolAccepted<=poolVisited && affectedAccepted<=affected,
						"Local response counts exceed their denominators");
					require(accepted==Math.addExact(poolAccepted,affectedAccepted),
						"Local response accepted count disagrees with pool and affected counts");
					require(netGain==Math.subtractExact(accepted,affected),
						"Local response net gain disagrees with accepted minus affected");
					require(affectedLost==Math.subtractExact(affected,affectedAccepted),
						"Local response affected loss disagrees with affected minus accepted");
					require(candidates==Math.addExact(Math.addExact(alignments,lowHits),hydro),
						"Local response candidate accounting disagrees");
					require(accepted<=alignments,"Local response accepted count exceeds alignments");
					require(candidates<=Math.multiplyExact(split ? 2 : 1,Math.addExact(poolVisited,affected)),
						"Local response candidate count exceeds visited tails times replacement models");
					p.response=(double)netGain;
				}else{
					for(int i=12; i<=27; i++){require(row.termEquals("NA",i),"Rejected local response retained fitted fields");}
				}
				rows.add(p);
			}
		}finally{close(reader,localResults);}
		require(!rows.isEmpty(),"Empty local response table");
		return rows;
	}

	/** Reads the explicit native export without converting audited requests to feasible ones. */
	private ArrayList<Probe> readFeatureTable(){
		final ArrayList<Probe> rows=new ArrayList<Probe>();
		final IntHashMap ids=new IntHashMap();
		final ByteFile reader=ByteFile.makeByteFile(featureTable,false);
		try{
			header(reader.nextLine(),FEATURE_HEADER,"local feature");
			final LineParser1 row=new LineParser1('\t');
			for(byte[] line=reader.nextLine();line!=null;line=reader.nextLine()){
				row.set(line); require(row.terms()==32,"Local feature rows require exactly 32 columns");
				require(row.termEquals("local_features_v1",0),"Unsupported local feature kind");
				final int id=identifier(row,1,false);
				require(ids.get(id)<0,"Duplicate candidate_index in feature table"); ids.put(id,rows.size());
				final long seed=decimal(row,2);
				final boolean split=row.termEquals("fission",3);
				require(split || row.termEquals("fusion",3),"Unknown operation in feature table");
				final int a=identifier(row,4,false),b=identifier(row,5,true);
				require(split ? b==-1 : b>=0 && b!=a,"Feature endpoints disagree with operation");
				final String status;
				if(row.termEquals("feasible_control",6)){status="feasible_control";}
				else if(row.termEquals("audited",6)){status="audited";}
				else{require(row.termEquals("unavailable",6),"Unknown feature request status");status="unavailable";}
				final boolean available=row.termEquals("available",7);
				require(available || row.termEquals("unavailable",7),"Unknown feature availability status");
				require(!available || !status.equals("unavailable"),"Unavailable request cannot have available fitted features");
				require(row.setBounds(8)>0,"Missing feature reason");
				final long members=row.termEquals("NA",9) ? 0 : nonnegative(row,9);
				nonnegative(row,31);
				final Probe p=new Probe(id,a,b,split,status.equals("feasible_control"),members);
				p.status=status; p.sampleSeed=Long.toString(seed); p.fitEligible=available; p.seen=true;
				p.responseAvailable=!row.termEquals("NA",20);
				if(available){
					require(!row.termEquals("NA",18) && !row.termEquals("NA",19) && p.responseAvailable,
						"Available features require both predictors and the measured response");
					p.cohesion=finite(row,18); p.separation=finite(row,19); p.response=finite(row,20);
				}else{
					// Explicitly unavailable rows survive in residuals and census, never as zero-valued observations.
					if(!row.termEquals("NA",18)){finite(row,18);}
					if(!row.termEquals("NA",19)){finite(row,19);}
					if(p.responseAvailable){p.response=finite(row,20);}
				}
				rows.add(p);
			}
		}finally{close(reader,featureTable);}
		require(!rows.isEmpty(),"Empty feature table");
		return rows;
	}

	/*--------------------------------------------------------------*/
	/*----------------       In-memory live fit     ----------------*/
	/*--------------------------------------------------------------*/

	/** Fits an intercept plus exactly three already-weighted features, without I/O.
	 * Every supplied row is an actual available observation; no missing-value or
	 * operation-type term is inferred. Inputs are read-only for the duration of this
	 * call and are never retained. Long measured gains are converted to doubles for
	 * fitting only; the selector/publisher retain their exact integer accounting.
	 * Uses the existing column-L2-normalized normal equations and pivoted solver;
	 * rank-deficient or numerically singular designs fail at its existing tolerance.
	 */
	public static LiveFit fitLive(final double[] cohesion,final double[] separation,
			final double[] sizePressure,final long[] measuredGain){
		require(cohesion!=null && separation!=null && sizePressure!=null && measuredGain!=null &&
			cohesion.length==measuredGain.length &&
			separation.length==measuredGain.length && sizePressure.length==measuredGain.length,
			"Live regression requires at least four aligned observations and all three feature arrays");
		final int n=measuredGain.length;
		final double[][] columns={cohesion,separation,sizePressure}, xtx=new double[4][4];
		final double[] scales=new double[4], xty=new double[4], x=new double[4];
		scales[0]=Math.sqrt(n);
		for(int term=0;term<columns.length;term++){
			for(double value:columns[term]){
				require(Double.isFinite(value),"Live regression features must be finite, never missing-value sentinels");
				scales[term+1]=Math.hypot(scales[term+1],value);
			}
			require(Double.isFinite(scales[term+1]),
				"Live regression feature column has zero or nonfinite norm");
		}
		if(n<4){throw new FitUnavailableException("insufficient_observations",
			"Live regression requires at least four aligned observations and all three feature arrays");}
		for(int term=1;term<scales.length;term++){
			if(scales[term]==0){throw new FitUnavailableException("singular_design","Live regression feature column has zero or nonfinite norm");}
		}
		double sumY=0;
		boolean constant=true;
		for(int row=0;row<n;row++){
			x[0]=1/scales[0];
			for(int term=0;term<columns.length;term++){x[term+1]=columns[term][row]/scales[term+1];}
			final double y=measuredGain[row];
			sumY+=y; constant&=measuredGain[row]==measuredGain[0];
			for(int i=0;i<x.length;i++){
				xty[i]+=x[i]*y;
				for(int j=0;j<x.length;j++){xtx[i][j]+=x[i]*x[j];}
			}
		}
		final double[] coefficients=solve(xtx,xty,true);
		for(int i=0;i<coefficients.length;i++){
			coefficients[i]/=scales[i];
			require(Double.isFinite(coefficients[i]),"Rescaling produced a nonfinite live regression coefficient");
		}
		// Rank was checked above even for constant responses. Their exact optimum
		// is the constant intercept, avoiding a spurious nonzero SSE from roundoff.
		if(constant){Arrays.fill(coefficients,0); coefficients[0]=measuredGain[0];}
		final double meanY=sumY/n;
		double sse=0,sst=0;
		for(int row=0;row<n;row++){
			final double predicted=predictLive(coefficients,cohesion[row],separation[row],sizePressure[row]);
			final double residual=measuredGain[row]-predicted, deviation=measuredGain[row]-meanY;
			sse+=residual*residual; sst+=deviation*deviation;
		}
		final double rmse=Math.sqrt(sse/n), r2=constant ? 1 : 1-sse/sst;
		require(Double.isFinite(meanY) && Double.isFinite(sse) && Double.isFinite(sst) &&
			Double.isFinite(rmse) && Double.isFinite(r2),"Live regression produced nonfinite residual metrics");
		return new LiveFit(n,coefficients,sse,rmse,r2);
	}

	/** A finite, well-formed live design lacks enough independent observations.
	 * Only diagnostic callers may recover from this specific exception. Malformed
	 * inputs, nonfinite arithmetic and coefficient overflow remain ordinary failures.
	 */
	public static final class FitUnavailableException extends IllegalArgumentException{
		private FitUnavailableException(final String reason_,final String message){
			super(message);reason=reason_;
		}
		public String reason(){return reason;}
		private final String reason;
		private static final long serialVersionUID=1L;
	}

	/** Immutable raw-unit coefficients and in-sample residual metrics, not a held-out accuracy claim. */
	public static final class LiveFit{
		private LiveFit(final int observations_,final double[] coefficients_,final double sse_,final double rmse_,final double r2_){
			assert(observations_>=4 && coefficients_.length==4) : "A live fit has an intercept and exactly three independently fitted features";
			observations=observations_; coefficients=coefficients_.clone(); sse=sse_; rmse=rmse_; r2=r2_;
		}
		public int observations(){return observations;}
		public double intercept(){return coefficients[0];}
		public double cohesionCoefficient(){return coefficients[1];}
		public double separationCoefficient(){return coefficients[2];}
		public double sizeCoefficient(){return coefficients[3];}
		public double sse(){return sse;}
		public double rmse(){return rmse;}
		/** Constant measured targets use their exact intercept solution, SSE zero and R-squared one. */
		public double rSquared(){return r2;}
		public double predict(final double q,final double s,final double p){return predictLive(coefficients,q,s,p);}
		private final int observations;
		private final double[] coefficients;
		private final double sse,rmse,r2;
	}

	/** Keeps fitted-row diagnostics and later caller predictions on the identical finite arithmetic path. */
	private static double predictLive(final double[] coefficients,final double q,final double s,final double p){
		require(Double.isFinite(q) && Double.isFinite(s) && Double.isFinite(p),"Live prediction requires three finite features");
		final double predicted=coefficients[0]+coefficients[1]*q+coefficients[2]*s+coefficients[3]*p;
		require(Double.isFinite(predicted),"Live regression prediction overflowed");
		return predicted;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Existing file fit      ----------------*/
	/*--------------------------------------------------------------*/

	private static Fit fit(final ArrayList<Probe> rows,final boolean fitSize,final boolean fitSlot){
		final int terms=2+(fitSize ? 1 : 0)+(fitSlot ? 1 : 0);
		require(rows.size()>=terms,"Ordinary least squares needs at least as many observations as fitted terms");
		final double[][] xtx=new double[terms][terms];
		final double[] xty=new double[terms];
		final double[] x=new double[terms], scales=new double[terms];
		double sumY=0;
		for(Probe p:rows){
			p.vector(x,fitSize,fitSlot);
			for(int i=0; i<terms; i++){scales[i]=Math.hypot(scales[i],x[i]);}
		}
		for(int i=0; i<terms; i++){
			require(scales[i]>0 && Double.isFinite(scales[i]),"Probe feature column is empty or nonfinite");
		}
		for(Probe p:rows){
			p.vector(x,fitSize,fitSlot);
			final double y=p.response;
			require(Double.isFinite(y),"Regression target overflowed");
			sumY+=p.response;
			for(int i=0; i<terms; i++){x[i]/=scales[i];}
			for(int i=0; i<terms; i++){
				xty[i]+=x[i]*y;
				for(int j=0; j<terms; j++){xtx[i][j]+=x[i]*x[j];}
			}
		}
		final double[] coefficients=solve(xtx,xty);
		for(int i=0; i<coefficients.length; i++){coefficients[i]/=scales[i];}
		final double meanY=sumY/rows.size();
		double sse=0,sst=0;
		for(Probe p:rows){
			p.vector(x,fitSize,fitSlot);
			final double predicted=dot(coefficients,x);
			final double residual=p.response-predicted;
			require(Double.isFinite(residual),"Regression residual overflowed");
			p.predicted=predicted; p.residual=residual;
			sse+=residual*residual;
			final double dev=p.response-meanY; sst+=dev*dev;
		}
		final double rmse=Math.sqrt(sse/rows.size());
		final double r2=sst==0 ? (sse==0 ? 1 : 0) : 1-sse/sst;
		require(Double.isFinite(sse) && Double.isFinite(rmse) && Double.isFinite(r2),
			"Probe regression produced nonfinite residual metrics");
		return new Fit(coefficients,sse,rmse,r2);
	}

	private static double[] solve(final double[][] matrix,final double[] rhs){
		return solve(matrix,rhs,false);
	}

	/** The live diagnostic path distinguishes finite rank deficiency from corrupt arithmetic. */
	private static double[] solve(final double[][] matrix,final double[] rhs,final boolean live){
		final int n=rhs.length;
		for(int pivot=0; pivot<n; pivot++){
			int best=pivot;
			for(int row=pivot+1; row<n; row++){
				if(Math.abs(matrix[row][pivot])>Math.abs(matrix[best][pivot])){best=row;}
			}
			if(live){
				require(Double.isFinite(matrix[best][pivot]),"Live regression pivot is nonfinite");
				if(Math.abs(matrix[best][pivot])<=1e-12){throw new FitUnavailableException("singular_design","Probe feature matrix is singular");}
			}else{require(Math.abs(matrix[best][pivot])>1e-12,"Probe feature matrix is singular");}
			if(best!=pivot){
				final double[] swap=matrix[pivot]; matrix[pivot]=matrix[best]; matrix[best]=swap;
				final double tmp=rhs[pivot]; rhs[pivot]=rhs[best]; rhs[best]=tmp;
			}
			final double scale=matrix[pivot][pivot];
			for(int col=pivot; col<n; col++){matrix[pivot][col]/=scale;}
			rhs[pivot]/=scale;
			for(int row=0; row<n; row++){
				if(row==pivot){continue;}
				final double factor=matrix[row][pivot];
				if(factor==0){continue;}
				for(int col=pivot; col<n; col++){matrix[row][col]-=factor*matrix[pivot][col];}
				rhs[row]-=factor*rhs[pivot];
			}
		}
		for(double value:rhs){require(Double.isFinite(value),"Probe regression produced a nonfinite coefficient");}
		return rhs;
	}

	private void write(final Fit fit,final int observations,final int requested){
		final ByteStreamWriter writer=new ByteStreamWriter(out,overwrite,false,true); writer.start();
		try{
			final ByteBuilder row=new ByteBuilder();
			writer.println("metric\tvalue");
			write(row,writer,"observations",observations);
			if(featureTable!=null){
				write(row,writer,"requested",requested); write(row,writer,"fit",observations);
				write(row,writer,"unavailable",requested-observations);
			}
			write(row,writer,"predictors",fit.coefficients.length);
			write(row,writer,"fit_slot",fitSlot ? 1 : 0);
			write(row,writer,"cohesion_scale",fit.coefficients[0]);
			write(row,writer,"separation_scale",fit.coefficients[1]);
			write(row,writer,"size_scale",fitSize ? fit.coefficients[2] : 0);
			write(row,writer,"slot_scale",fitSlot ? fit.coefficients[fitSize ? 3 : 2] : 1);
			write(row,writer,"sse",fit.sse);
			write(row,writer,"rmse",fit.rmse);
			write(row,writer,"r_squared",fit.r2);
		}finally{if(writer.poisonAndWait()){throw new IllegalStateException("Failed probe-regression output: "+out);}}
	}

	/** Preserves every request and its status; only fitted rows have predictions. */
	private void writeResiduals(final ArrayList<Probe> probes){
		final ByteStreamWriter writer=new ByteStreamWriter(outResiduals,overwrite,false,true); writer.start();
		try{
			writer.println(RESIDUAL_HEADER);
			final ByteBuilder row=new ByteBuilder();
			for(Probe p:probes){
				row.clear().append(p.id).tab().append(p.sampleSeed).tab().append(p.split ? "fission" : "fusion")
					.tab().append(p.a).tab().append(p.b).tab().append(p.status).tab();
				if(p.fitEligible){row.appendSlow(p.response).tab().appendSlow(p.predicted).tab().appendSlow(p.residual);}
				else{
					if(p.responseAvailable){row.appendSlow(p.response);}else{row.append("NA");}
					row.append("\tNA\tNA");
				}
				writer.println(row);
			}
		}finally{if(writer.poisonAndWait()){throw new IllegalStateException("Failed residual output: "+outResiduals);}}
	}

	private static void register(final ArrayList<Path> paths,final String name,final boolean input) throws IOException{
		require(name!=null && !name.isEmpty() && name.indexOf('\t')<0 && name.indexOf('\n')<0 && name.indexOf('\r')<0,
			"Regression paths must be nonempty unescaped TSV fields");
		final Path path=Paths.get(name).toAbsolutePath().normalize();
		require(path.getFileName()!=null,"Regression paths must name disk files");
		final Path real;
		if(input || Files.exists(path,LinkOption.NOFOLLOW_LINKS)){
			real=path.toRealPath(); require(Files.isRegularFile(real),"Regression input/output must be a regular file");
		}else{
			Path ancestor=path.getParent();
			while(ancestor!=null && !Files.exists(ancestor,LinkOption.NOFOLLOW_LINKS)){ancestor=ancestor.getParent();}
			require(ancestor!=null,"Path has no existing ancestor");
			real=ancestor.toRealPath().resolve(ancestor.relativize(path)).normalize();
		}
		for(Path old:paths){
			require(!old.equals(real) && !(Files.exists(old) && Files.exists(real) && Files.isSameFile(old,real)),
				"Probe-regression source/output paths alias");
		}
		paths.add(real);
	}
	private static long decimal(final LineParser1 row,final int term){
		require(row.setBounds(term)>0,"Missing integer field");
		final byte[] line=row.line(); int at=row.a(); final int end=row.b();
		final boolean negative=line[at]=='-'; if(negative){at++;}
		require(at<end,"Missing integer digits");
		final long limit=negative ? Long.MIN_VALUE : -Long.MAX_VALUE; long value=0;
		for(; at<end; at++){
			final int digit=line[at]-'0';
			require(digit>=0 && digit<=9 && value>=limit/10,"Invalid or overflowing integer");
			value*=10; require(value>=limit+digit,"Overflowing integer"); value-=digit;
		}
		return negative ? value : -value;
	}
	private static long nonnegative(final LineParser1 row,final int term){
		final long value=decimal(row,term); require(value>=0,"Negative count field"); return value;
	}
	private static int identifier(final LineParser1 row,final int term,final boolean absent){
		final long value=decimal(row,term);
		require(value>=(absent ? -1 : 0) && value<Integer.MAX_VALUE,"Identifier out of range"); return (int)value;
	}
	private static double finite(final LineParser1 row,final int term){
		require(row.setBounds(term)>0,"Missing numeric field");
		final double value=Parse.parseDoubleSlow(row.line(),row.a(),row.b());
		require(Double.isFinite(value),"Nonfinite numeric field"); return value;
	}
	private static void header(final byte[] line,final String expected,final String description){
		require(line!=null && Arrays.equals(line,expected.getBytes(java.nio.charset.StandardCharsets.US_ASCII)),
			"Unexpected "+description+" header");
	}
	private static void close(final ByteFile reader,final String path){
		if(reader.close()){throw new IllegalStateException("Failed probe-regression input: "+path);}
	}
	private static void write(final ByteBuilder row,final ByteStreamWriter writer,final String key,final int value){
		row.clear().append(key).tab().append(value); writer.println(row);
	}
	private static void write(final ByteBuilder row,final ByteStreamWriter writer,final String key,final double value){
		row.clear().append(key).tab().appendSlow(value); writer.println(row);
	}
	private static double dot(final double[] a,final double[] b){
		double sum=0; for(int i=0; i<a.length; i++){sum+=a[i]*b[i];} return sum;
	}
	private static void require(final boolean valid,final String message){if(!valid){throw new IllegalArgumentException(message);}}

	private static final class Probe{
		Probe(final int id_,final int a_,final int b_,final boolean split_,
				final boolean feasible_,final long initial_){
			id=id_; a=a_; b=b_; split=split_; feasible=feasible_; initial=initial_;
			fitEligible=feasible_; responseAvailable=feasible_; status=feasible_ ? "feasible" : "infeasible";
		}
		void vector(final double[] x,final boolean fitSize,final boolean fitSlot){
			x[0]=cohesion; x[1]=separation;
			int next=2; if(fitSize){x[next++]=size;} if(fitSlot){x[next]=slot;}
		}
		final int id,a,b;
		final long initial;
		final boolean split,feasible;
		boolean seen,fitEligible,responseAvailable;
		String sampleSeed="NA",status;
		double cohesion,separation,size,slot,response,predicted,residual;
	}
	private static final class Fit{
		Fit(final double[] coefficients_,final double sse_,final double rmse_,final double r2_){
			coefficients=coefficients_; sse=sse_; rmse=rmse_; r2=r2_;
		}
		final double[] coefficients;
		final double sse,rmse,r2;
	}

	public static final String RESIDUAL_HEADER="candidate_index\tsample_seed\toperation\tleft_id\tright_id\tstatus\tmeasured\tpredicted\tresidual";
	public static final String FEATURE_HEADER=
		"report_kind\tcandidate_index\tsample_seed\toperation\tleft_id\tright_id\tstatus\tfeature_status\treason"+
		"\tparent_members\tleft_members\tright_members\tparent_q\tparent_n\tleft_q\tleft_n\tright_q\tright_n"+
		"\tweighted_cohesion_delta\tweighted_separation_delta\tnet_gain"+
		"\tleft_pair_available\tleft_pair_identity\tleft_pair_score_ratio\tleft_pair_shared_fraction\tleft_pair_quality"+
		"\tright_pair_available\tright_pair_identity\tright_pair_score_ratio\tright_pair_shared_fraction\tright_pair_quality\telapsed_ns";
	private String results,report,localResults,featureTable,out,outResiduals;
	private boolean overwrite,fitSlot;
	private boolean fitSize=true;
	private final PrintStream outstream;
}
