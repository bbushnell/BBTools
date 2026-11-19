package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Locale;
import java.util.TreeSet;

import ddl.DDLPentamerBlacklist;
import ddl.DDLFamilySearch.StructuralInputs;
import ddl.ProteinKmerGeometry;
import fileIO.ByteFile;
import parse.LineParser1;
import parse.Parse;
import shared.Resources;
import structures.ByteBuilder;

/** Resolves one classifier recipe without changing its scoring policy.
 * @author Brian, Keqing
 */
public final class FamilyScoringPackage{

	/** Bind a recipe file or a directory containing classifier.config; validate component pins. */
	public FamilyScoringPackage(Path source) throws Exception{
		this(source,FamilyCalibration.VerificationMode.SHA80);
	}

	/** Explicit structural loading trusts declared byte identities while retaining schema,
	 * roster and observable file-stability checks, including explicit filtered packages.
	 * Construction, loading and verification share single-threaded tracking state;
	 * loaded classifiers still provide their ordinary independent workers. */
	public FamilyScoringPackage(Path source,FamilyCalibration.VerificationMode verificationMode_) throws Exception{
		require(source!=null && verificationMode_!=null,"Package binding needs a path and explicit verification policy");
		verificationMode=verificationMode_;
		structuralInputs=verificationMode==FamilyCalibration.VerificationMode.STRUCTURAL ?
			new StructuralInputs() : null;
		final CompactSource compact=compactFamily(source,structuralInputs);
		compactFamily=compact==null ? null : compact.family;
		compactReference=compact==null ? null : compact.reference;
		if(compact!=null){
			require(structuralInputs!=null,
				"Compact package resources carry structural identity; use explicit structural loading");
			compactRoot=compact.family.getParent().getParent().toAbsolutePath().normalize();
			path=compact.family.toAbsolutePath().normalize();
			if(Files.notExists(path,LinkOption.NOFOLLOW_LINKS)){throw missingCompactFamily(path);}
			structuralInputs.track(path);
			fingerprint="unmeasured";
			compactHeader=readCompactHeader(path);
			bindCompactPackage();
			gcInputs=null; memberRange=MemberLengthRange.from(options);
			matrix=resolveMatrix(options);
			structuralInputs.verify();
			return;
		}
		compactRoot=null; compactHeader=null;
		path=(Files.isDirectory(source) ? source.resolve("classifier.config") : source).toAbsolutePath().normalize();
		if(Files.notExists(path,LinkOption.NOFOLLOW_LINKS)){throw missingLegacyComponent(path,"recipe");}
		require(Files.isRegularFile(path),"Missing classifier package recipe: "+path);
		if(structuralInputs!=null){structuralInputs.track(path);fingerprint="unmeasured";}
		else{fingerprint=DigestSuffix.file(path.toString());}
		final ByteFile input=ByteFile.makeByteFile(path.toString(),false);
		try{
			for(byte[] line=input.nextLine();line!=null;line=input.nextLine()){
				final String text=new String(line,StandardCharsets.UTF_8).trim();
				if(text.isEmpty() || text.startsWith("#")){continue;}
				final String[] pair=Parse.splitOnFirst(text,'=');
				require(pair.length==2,"Package recipe requires key=value: "+text);
				final String key=pair[0].trim().toLowerCase(Locale.ROOT),value=pair[1].trim();
				require(!key.isEmpty() && !value.isEmpty() && options.put(key,value)==null,
					"Empty or duplicate package option: "+key);
			}
		}finally{require(!input.close(),"Failed reading classifier package: "+path);}
		final String declaredVerification=options.remove("verification");
		require(declaredVerification==null || structuralInputs!=null && declaredVerification.equals("structural"),
			"An inherited-identity package requires explicit structural loading; default loading cannot certify it");
		normalizeStructuralLabels(options);
		FamilyBlacklistBinding.requireOfficialOptions(options);
		FamilyBlacklistBinding.validateOptions(options,verificationMode);
		require(structuralInputs==null || !options.containsKey("blacklist") ||
			FamilyBlacklistBinding.structuralFileGeometry(options),
			"Structural file-backed filtering requires a supported explicit full-word geometry");
		if(options.containsKey("mode")){
			require("search".equals(options.remove("mode")),"Evaluation requires an ordinary search package");
		}
		for(String key:REQUIRED){require(options.containsKey(key),"Package must explicitly specify "+key);}
		validateRankingOptions(options,verificationMode);
		gcInputs=options.containsKey("familygc") ?
			(structuralInputs==null ? new StructuralInputs() : structuralInputs) : null;
		memberRange=MemberLengthRange.from(options);
		for(String key:PATHS){
			if(!options.containsKey(key)){
				require(!options.containsKey(key+"sha80"),"Package fingerprint has no path: "+key);
				continue;
			}
			final Path value=path.getParent().resolve(options.get(key)).normalize();
			if(Files.notExists(value,LinkOption.NOFOLLOW_LINKS)){throw missingLegacyComponent(value,key);}
			require(Files.isRegularFile(value),"Missing package component "+key+": "+value);
			options.put(key,value.toString());
			if(gcInputs!=null && (key.equals("familygc") || key.equals("familygcmap"))){gcInputs.track(value);}
			final String expected=options.get(key+"sha80");
			if(structuralInputs!=null){
				if(expected!=null && !"unmeasured".equals(expected)){
					DigestSuffix.requireSuffix(expected,"Inherited package component "+key);
				}
				structuralInputs.track(value);
			}else{
				final String actual=DigestSuffix.file(value.toString());
				require(expected==null || expected.equals(actual),"Package component fingerprint differs: "+key);
				options.put(key+"sha80",actual);
			}
		}
		matrix=resolveMatrix(options);
		if(structuralInputs!=null){structuralInputs.verify();}
		if(gcInputs!=null && gcInputs!=structuralInputs){gcInputs.verify();}
	}

	/** Structural mode never establishes or recomputes payload byte identity. */
	public FamilyCalibration.VerificationMode verificationMode(){return verificationMode;}

	/** Detached resolved metadata for resource discovery/output protection. Loading
	 * must still use this package; editing this copy never changes its policy. */
	public HashMap<String,String> resolvedOptions(){return new HashMap<String,String>(options);}

	/** Immutable compact metadata, or null for a legacy recipe. Reading it does not load families. */
	public FamilyCompactCodec.Header compactHeader(){return compactHeader;}
	/** Detached roster of the compact family file and its three resolved shared resources.
	 * Callers must protect these files from output aliases and retain verify() checks. */
	public Path[] compactInputs(){
		require(compactFamily!=null && compactInputs!=null,"Compact input discovery requires a compact package");
		return compactInputs.clone();
	}
	/** Reloadable reference pointer used to reach a compact family, or null for direct/legacy sources. */
	public Path compactReference(){return compactReference;}

	/** Compact evaluation outputs must not alias or nest under the immutable package. */
	public void protectCompactOutput(Path output) throws Exception{
		if(compactFamily==null){return;}
		require(output!=null,"Compact output protection requires a destination");
		require(!Files.exists(output,LinkOption.NOFOLLOW_LINKS),
			"Compact evaluation output directory already exists: "+output);
		final Path rootReal=compactRoot.toRealPath();
		Path existing=output.toAbsolutePath().normalize().getParent();
		while(existing!=null && !Files.exists(existing,LinkOption.NOFOLLOW_LINKS)){
			existing=existing.getParent();
		}
		require(existing!=null,"Compact evaluation output needs an existing parent");
		final Path parentReal=existing.toRealPath();
		require(!parentReal.startsWith(rootReal),
			"Compact evaluation output cannot write inside the compact package root: "+output);
	}

	/** Versioned pointer to the original compact family file; effective run policy is recorded separately. */
	ByteBuilder compactReferenceRecipe(){
		require(compactFamily!=null,"Only compact packages emit compact references");
		return new ByteBuilder(COMPACT_REFERENCE_MAGIC).tab().append(COMPACT_REFERENCE_VERSION).nl()
			.append("view=").append(required(compactHeader,"view")).nl()
			.append("family=").append(path.toString()).nl();
	}

	/** Load the same immutable classifier used by FamilyEval, without applying a margin override. */
	public FamilyTailCascade.LoadedClassifier load(int threads) throws Exception{
		assert(threads>0) : "Each classifier worker needs a positive configured thread budget";
		return loadSelected(options,threads);
	}

	/** Explicit shortlist comparison; clone the recipe before loading so the source policy stays immutable. */
	public FamilyTailCascade.LoadedClassifier loadWithTop(int threads,int top) throws Exception{
		require(top>0 && threads>0,"An explicit shortlist cap and worker count must be positive");
		return loadWithOverrides(threads,top,-1);
	}
	/** Independent fixed-floor/cap experiments copy options before loading; neither edits the parent recipe. */
	public FamilyTailCascade.LoadedClassifier loadWithOverrides(int threads,int top,int ddlReject) throws Exception{
		require(threads>0 && top>=0 && ddlReject>=-1,
			"Overrides need positive workers, optional positive cap and optional nonnegative raw-hit floor");
		final HashMap<String,String> selected=new HashMap<String,String>(options);
		if(top>0){selected.put("top",Integer.toString(top));}
		if(ddlReject>=0){selected.put("ddlreject",Integer.toString(ddlReject));}
		return loadSelected(selected,threads);
	}

	/** Explicit evaluator-only policy; preserve the package and every non-GC scoring option.
	 * FamilyEval records this override separately from its original-package reference.
	 */
	FamilyTailCascade.LoadedClassifier loadForEvaluationWithoutGc(int threads,int top,int ddlReject)throws Exception{
		require(threads>0 && top>=0 && ddlReject>=-1,
			"Evaluation overrides need positive workers, optional positive cap and optional nonnegative raw-hit floor");
		final HashMap<String,String> selected=new HashMap<String,String>(options);
		if(top>0){selected.put("top",Integer.toString(top));}
		if(ddlReject>=0){selected.put("ddlreject",Integer.toString(ddlReject));}
		for(String key:GC_RANKING){selected.remove(key);}
		selected.remove("gcfilter");selected.remove("gc_z");
		final FamilyTailCascade.LoadedClassifier result=loadSelected(selected,threads);
		require(!result.usesFamilyGcPenalty(),"An explicit GC-off evaluation must not retain GC ranking or filtering");
		return result;
	}

	/** Preserve default loading exactly; propagate the explicit structural context through every file reader. */
	private FamilyTailCascade.LoadedClassifier loadSelected(HashMap<String,String> selected,int threads)throws Exception{
		if(compactFamily!=null){
			structuralInputs.verify();
			final FamilyCompactLoader.Loaded loaded=FamilyCompactLoader.load(compactFamily.toString(),threads,1000000);
			final ProteinKmerGeometry geometry=ProteinKmerGeometry.make(Integer.parseInt(selected.get("k")),
				selected.get("aabits"));
			final DDLPentamerBlacklist blacklist=DDLPentamerBlacklist.loadReducedStructural(
				Paths.get(selected.get("blacklist")),"unmeasured",
				Integer.parseInt(selected.get("blacklistwords")),geometry,selected.get("blacklistprovenance"));
			final FamilyTailCascade.LoadedClassifier result=
				FamilyTailCascade.loadCompactCallable(selected,threads,loaded,blacklist);
			structuralInputs.verify(); return result;
		}
		if(gcInputs!=null){gcInputs.verify();}
		if(structuralInputs==null){
			return bindRanking(FamilyTailCascade.loadCallable(baseOptions(selected),threads),selected,verificationMode,gcInputs);
		}
		structuralInputs.verify();
		final FamilyTailCascade.LoadedClassifier base=FamilyTailCascade.loadCallableStructural(
			baseOptions(selected,verificationMode),threads,structuralInputs);
		final FamilyTailCascade.LoadedClassifier result=bindRanking(base,selected,verificationMode,gcInputs);
		structuralInputs.verify();return result;
	}

	/** Structural package recipes may omit component byte pins; make the
	 * unmeasured identity explicit before downstream loaders require a label for
	 * every supplied file. */
	private void normalizeStructuralLabels(HashMap<String,String> source){
		if(structuralInputs==null){return;}
		for(String key:PATHS){
			if(source.containsKey(key) && !source.containsKey(key+"sha80")){
				source.put(key+"sha80","unmeasured");
			}
		}
	}

	/** Removes only the wrapper binding; alignment mode/span remain loader-validated options. */
	static HashMap<String,String> baseOptions(HashMap<String,String> source){
		return baseOptions(source,FamilyCalibration.VerificationMode.SHA80);
	}
	/** Resolve immutable score units from exactly one resource or portable binding.
	 * An absent selector preserves historical BLOSUM62; a named custom matrix cannot fall back.
	 */
	public static ProteinSubstitutionMatrix resolveMatrix(HashMap<String,String> source){
		require(source!=null,"Matrix selection requires explicit options");
		final String name=source.get("matrix"),file=source.get("matrixfile"),binding=source.get("matrixbinding");
		require(file==null || binding==null,"Choose matrixfile or matrixbinding, not both");
		require(file!=null || !source.containsKey("matrixfilesha80"),"Matrix fingerprint requires a resource path");
		if(name==null){
			require(file==null && binding==null,"Matrix resource/binding requires an explicit matrix selector");
			return ProteinSubstitutionMatrix.BLOSUM62;
		}
		require(ProteinSubstitutionMatrix.supportedSelector(name),"Unknown matrix selector: "+name);
		if(binding!=null){
			final ProteinSubstitutionMatrix result=ProteinSubstitutionMatrix.parseBinding(binding);
			require(result.name.equals(name),"Matrix selector and embedded scores disagree");return result;
		}
		return ProteinSubstitutionMatrix.load(name,file);
	}
	/** New preparation selects the current default and records portable scores explicitly.
	 * Call after preparation has pinned external matrix inputs. Historical deserialization
	 * deliberately continues to call resolveMatrix, whose absent selector means BLOSUM62.
	 */
	static ProteinSubstitutionMatrix resolveNewWorkMatrix(HashMap<String,String> source){
		require(source!=null,"New matrix selection requires explicit options");
		if(!source.containsKey("matrix")){
			require(!source.containsKey("matrixfile") && !source.containsKey("matrixfilesha80") &&
				!source.containsKey("matrixbinding"),"Matrix resource/binding requires an explicit matrix selector");
			source.put("matrix",ProteinSubstitutionMatrix.DEFAULT_NAME);
		}
		final ProteinSubstitutionMatrix matrix=resolveMatrix(source);
		source.remove("matrixfile");source.remove("matrixfilesha80");
		source.put("matrixbinding",matrix.binding());
		return matrix;
	}
	/** Explicit structural loading validates labels before removing the wrapper-only ranking fields. */
	static HashMap<String,String> baseOptions(HashMap<String,String> source,FamilyCalibration.VerificationMode mode){
		validateRankingOptions(source,mode);
		final HashMap<String,String> result=new HashMap<String,String>(source);
		for(String key:RANKING){result.remove(key);}
		for(String key:GC_RANKING){result.remove(key);}
		return result;
	}

	/** Compact Q76 packages use one families.tsv per view plus shared package-root resources. */
	private void bindCompactPackage() throws Exception{
		require(Files.isRegularFile(path,LinkOption.NOFOLLOW_LINKS),
			"Missing compact classifier family file: "+path);
		final String view=path.getParent().getFileName().toString();
		require(("short".equals(view) || "tails".equals(view)) && view.equals(compactHeader.get("view")),
			"Compact family view must match its short/ or tails/ directory");
		final Path rootReal=compactRoot.toRealPath();
		final Path alphabet=resolveCompactResource(rootReal,"alphabet");
		final Path blacklist=resolveCompactResource(rootReal,"blacklist");
		final Path matrixResource=resolveCompactResource(rootReal,"matrixresource");
		structuralInputs.track(alphabet);
		structuralInputs.track(blacklist);
		structuralInputs.track(matrixResource);
		compactInputs=new Path[]{path,alphabet,blacklist,matrixResource};
		options.put("aabits",readCompactAlphabet(alphabet));
		copyCompactKeys("top","k","seed","width","ddlreject","ddlaccept","cosine","rankingfeatures",
			"scalarscale","scalargates","scalarz","identityacceptfloor","blosumacceptfloor",
			"partitions","overlap","blacklistk","blacklistmode","blacklistwords","blacklistprovenance",
			"indexearlyreject","alignmentmode","alignmentminspan","alignmentminpaired",
			"heuristic","minqueryaa","maxqueryaa","familygclambda","gcfilter","gc_z",
			"minmemberaa","maxmemberaa");
		options.put("positionmode","fixed");
		options.put("blacklist",blacklist.toString());
		options.put("matrix",required(compactHeader,"matrix"));
		final ProteinSubstitutionMatrix loaded=ProteinSubstitutionMatrix.load(options.get("matrix"),
			matrixResource.toString());
		require(loaded.binding().equals(required(compactHeader,"matrixbinding")),
			"Compact package matrix resource differs from the embedded binding");
		options.put("matrixbinding",loaded.binding());
	}

	/** Shared locators are package-root relative and may not escape by syntax or symlink. */
	private Path resolveCompactResource(Path rootReal,String key) throws Exception{
		final String value=required(compactHeader,key);
		final Path relative=Paths.get(value);
		require(!relative.isAbsolute() && value.charAt(0)!='?',
			"Compact "+key+" must be relative to the package root");
		final Path resolved=compactRoot.resolve(relative).normalize();
		require(resolved.startsWith(compactRoot),"Compact "+key+" escapes the package root: "+value);
		if(Files.notExists(resolved,LinkOption.NOFOLLOW_LINKS)){
			// A missing child must not hide an existing parent symlink escaping the package.
			Path parent=resolved.getParent();
			while(parent!=null && Files.notExists(parent,LinkOption.NOFOLLOW_LINKS)){parent=parent.getParent();}
			require(parent!=null && Files.isDirectory(parent) && parent.toRealPath().startsWith(rootReal),
				"Compact "+key+" does not resolve within the package root: "+value);
			throw new Resources.MissingResourceException("\nERROR: Required compact BBGAP "+key+
				" resource not found: "+resolved+"\nRestore the small resources supplied with BBTools under "+
				compactRoot+", preserving the filenames declared by the family header.\n"+
				"Obtain the matching BBTools release from https://github.com/bbushnell/BBTools/releases/\n");
		}
		final Path real=resolved.toRealPath();
		require(real.startsWith(rootReal) && Files.isRegularFile(real),
			"Compact "+key+" does not resolve to an in-package regular file: "+value);
		return real;
	}

	private static FamilyCompactCodec.Header readCompactHeader(Path file) throws Exception{
		final ByteFile input=ByteFile.makeByteFile1(file.toString(),false);
		final byte[][] lines=new byte[FamilyCompactCodec.LINES_PER_RECORD][];
		try{
			for(int i=0;i<lines.length;i++){
				lines[i]=input.nextLine();
				require(lines[i]!=null,"Incomplete compact header: "+file);
			}
			return FamilyCompactCodec.parseHeader(lines,0,new LineParser1('\t'));
		}finally{require(!input.close(),"Failed reading compact header: "+file);}
	}

	private static String readCompactAlphabet(Path file) throws Exception{
		final ByteFile input=ByteFile.makeByteFile1(file.toString(),false);
		try{
			final byte[] first=input.nextLine(),second=input.nextLine();
			require(first!=null && new String(first,StandardCharsets.US_ASCII).equals("#bbgap_alphabet_config\t1") &&
				second!=null,"Compact alphabet resource is missing its schema or mapping");
			final String[] pair=Parse.splitOnFirst(new String(second,StandardCharsets.US_ASCII),'=');
			require(pair.length==2 && pair[0].equals("aabits") && input.nextLine()==null,
				"Compact alphabet resource must contain one aabits mapping");
			return pair[1];
		}finally{require(!input.close(),"Failed reading compact alphabet: "+file);}
	}

	private void copyCompactKeys(String... keys){
		for(String key:keys){
			final String value=compactHeader.get(key);
			if(value!=null){options.put(key.toLowerCase(Locale.ROOT),value);}
		}
	}

	private static CompactSource compactFamily(Path source,StructuralInputs tracker)throws Exception{
		final Path path=source.toAbsolutePath().normalize();
		if(isCompactFamily(path) && Files.exists(path,LinkOption.NOFOLLOW_LINKS)){
			require(Files.isRegularFile(path),"Compact family file must be regular: "+path);
		}
		if(!Files.isDirectory(source)){
			if(isCompactFamily(path)){return new CompactSource(path,null);}
			if(compactView(path) && Files.notExists(path,LinkOption.NOFOLLOW_LINKS)){
				throw missingCompactFamily(path.resolve("families.tsv.gz"));
			}
			final CompactSource reference=compactReference(path,tracker);
			if(reference!=null){return reference;}
			final String name=path.getFileName().toString();
			require(!(name.equals("families.tsv") || name.equals("families.tsv.gz")),
				"Compact family files must be under short/ or tails/: "+path);
			return null;
		}
		final Path compressed=path.resolve("families.tsv.gz");
		if(compactView(path) && Files.isRegularFile(compressed)){return new CompactSource(compressed,null);}
		final Path plain=path.resolve("families.tsv");
		if(compactView(path) && Files.isRegularFile(plain)){return new CompactSource(plain,null);}
		if(compactView(path)){
			require(!Files.exists(compressed,LinkOption.NOFOLLOW_LINKS) && !Files.exists(plain,LinkOption.NOFOLLOW_LINKS),
				"Compact family view contains a nonregular family entry: "+path);
		}
		if(compactView(path) && Files.notExists(compressed,LinkOption.NOFOLLOW_LINKS) &&
			Files.notExists(plain,LinkOption.NOFOLLOW_LINKS) &&
			Files.notExists(path.resolve("classifier.config"),LinkOption.NOFOLLOW_LINKS)){
			throw missingCompactFamily(compressed);
		}
		final boolean hasShort=Files.exists(path.resolve("short").resolve("families.tsv.gz")) ||
			Files.exists(path.resolve("short").resolve("families.tsv"));
		final boolean hasTails=Files.exists(path.resolve("tails").resolve("families.tsv.gz")) ||
			Files.exists(path.resolve("tails").resolve("families.tsv"));
		require(!hasShort && !hasTails,
			"Compact package roots are ambiguous; select short/families.tsv.gz or tails/families.tsv.gz");
		return null;
	}

	/** Names the archive for an absent compact view; existing malformed files still reach their parser. */
	private static Resources.MissingResourceException missingCompactFamily(Path file){
		assert(isCompactFamily(file)) : "Archive guidance requires a recognized short/tails family-file path";
		final String view=file.getParent().getFileName().toString();
		final String archive=view.equals("short") ? "bbgap_short_v1.tar" : "bbgap_tails_v1.tar";
		return new Resources.MissingResourceException("\nERROR: Required BBGAP "+view+
			" family resource not found: "+file+"\n"+Resources.bbgapDownloadInstructions(archive));
	}

	/** Legacy recipes are user-produced artifacts, not files supplied by the compact-view archives. */
	private static Resources.MissingResourceException missingLegacyComponent(Path file,String role){
		return new Resources.MissingResourceException("\nERROR: Required classifier package "+role+
			" not found: "+file+"\nRestore or regenerate the complete legacy classifier package with its original producer,\n"+
			"including classifier.config and every component named by that recipe.\n"+
			"For distributed BBGAP families, obtain the compact SHORT or TAILS archive from\n"+
			"https://github.com/bbushnell/BBTools/releases/ and select its short/families.tsv.gz\n"+
			"or tails/families.tsv.gz with in=. Compact archives do not recreate a custom legacy package.\n");
	}

	private static CompactSource compactReference(Path source,StructuralInputs tracker)throws Exception{
		if(!Files.isRegularFile(source,LinkOption.NOFOLLOW_LINKS)){return null;}
		if(tracker!=null){tracker.track(source);}
		final ByteFile input=ByteFile.makeByteFile1(source.toString(),false);
		try{
			final byte[] first=input.nextLine();
			if(first==null){return null;}
			final String magic=new String(first,StandardCharsets.US_ASCII);
			if(magic.startsWith(COMPACT_REFERENCE_MAGIC+"\t")){
				require(magic.equals(COMPACT_REFERENCE_MAGIC+"\t"+COMPACT_REFERENCE_VERSION),
					"Unsupported compact reference version: "+magic);
			}else{return null;}
			require(tracker!=null,"Compact reference pointers require explicit structural loading");
			tracker.track(source);
			final HashMap<String,String> fields=new HashMap<String,String>();
			for(byte[] line=input.nextLine();line!=null;line=input.nextLine()){
				final String text=new String(line,StandardCharsets.US_ASCII);
				final String[] pair=Parse.splitOnFirst(text,'=');
				require(pair.length==2 && pair[0].length()>0 && pair[1].length()>0 &&
					fields.put(pair[0],pair[1])==null,"Compact reference requires unique key=value rows");
			}
			require(fields.size()==2 && fields.containsKey("view") && fields.containsKey("family"),
				"Compact reference requires exactly view and family rows");
			final String view=fields.get("view");
			require(view.equals("short") || view.equals("tails"),"Compact reference has an unsupported view: "+view);
			final Path target=Paths.get(fields.get("family"));
			require(target.isAbsolute() && target.normalize().equals(target),
				"Compact reference family must be a normalized absolute path");
			require(isCompactFamily(target) && target.getParent().getFileName().toString().equals(view),
				"Compact reference must target its exact compact view family file");
			return new CompactSource(target,source);
		}finally{require(!input.close(),"Failed reading compact reference: "+source);}
	}

	/** Recognizes view filenames without treating a filesystem root as a named file. */
	private static boolean isCompactFamily(Path path){
		if(path==null || path.getFileName()==null){return false;}
		final String name=path.getFileName().toString();
		return (name.equals("families.tsv") || name.equals("families.tsv.gz")) &&
			compactView(path.getParent());
	}

	private static boolean compactView(Path dir){
		if(dir==null || dir.getFileName()==null){return false;}
		final String name=dir.getFileName().toString();
		return name.equals("short") || name.equals("tails");
	}

	private static String required(FamilyCompactCodec.Header header,String key){
		final String value=header.get(key);
		require(value!=null && value.length()>0,"Compact package header missing "+key);
		return value;
	}

	/** Adds the recipe's pinned common pivot to the already loaded index; never loads a second classifier. */
	static FamilyTailCascade.LoadedClassifier bindRanking(FamilyTailCascade.LoadedClassifier base,
			HashMap<String,String> source)throws Exception{
		return bindRanking(base,source,FamilyCalibration.VerificationMode.SHA80);
	}
	/** Resident preparation explicitly propagates input verification without changing ranking policy. */
	static FamilyTailCascade.LoadedClassifier bindRanking(FamilyTailCascade.LoadedClassifier base,
			HashMap<String,String> source,FamilyCalibration.VerificationMode verificationMode)throws Exception{
		return bindRanking(base,source,verificationMode,null);
	}
	/** GC callers retain this tracker through final verification; no implicit untracked resident load. */
	static FamilyTailCascade.LoadedClassifier bindRanking(FamilyTailCascade.LoadedClassifier base,
			HashMap<String,String> source,FamilyCalibration.VerificationMode verificationMode,
			StructuralInputs gcTracker)throws Exception{
		require(verificationMode!=null,"Ranking needs an explicit verification policy");
		validateRankingOptions(source,verificationMode);
		if(!source.containsKey("heuristic")){return base;}
		final FamilyTailCascade.LoadedClassifier ranked=base.withHeuristicRanking(new FamilyHeuristicScore(base,
			java.nio.file.Paths.get(source.get("lengthmedians")),
			source.get("lengthmedianssha80"),java.nio.file.Paths.get(source.get("familyhydro")),source.get("familyhydrosha80"),
			Integer.parseInt(source.get("minqueryaa")),Integer.parseInt(source.get("maxqueryaa")),
			source.get("heuristic"),verificationMode));
		if(!source.containsKey("familygc")){return ranked;}
		require(gcTracker!=null,"GC package loading requires an explicitly retained input-stability tracker");
		if(verificationMode==FamilyCalibration.VerificationMode.SHA80){verifyGcIdentities(source);}
		final FamilyGcPenalty penalty=FamilyGcInputs.load(
			Paths.get(source.get("familygc")),Paths.get(source.get("familygcmap")),
			Integer.parseInt(source.get("familygcfamilies")),Long.parseLong(source.get("familygcmembers")),
			Double.parseDouble(source.get("familygclambda")),gcTracker);
		if(verificationMode==FamilyCalibration.VerificationMode.SHA80){verifyGcIdentities(source);}
		return ranked.withPackagedFamilyGcPenalty(penalty);
	}
	/** Strict loading verifies the recipe's content identities on both sides of the GC read. */
	private static void verifyGcIdentities(HashMap<String,String> source){
		for(String key:new String[]{"familygc","familygcmap"}){
			require(source.get(key+"sha80").equals(DigestSuffix.file(source.get(key))),"GC resource fingerprint differs: "+key);
		}
	}

	/** Complete bindings are explicit; partially specified or silently ignored ranking inputs are errors. */
	static void validateRankingOptions(HashMap<String,String> source){
		validateRankingOptions(source,FamilyCalibration.VerificationMode.SHA80);
	}
	/** Structural callers may declare new metadata unmeasured; omitted policy remains strict SHA80. */
	static void validateRankingOptions(HashMap<String,String> source,FamilyCalibration.VerificationMode mode){
		require(source!=null && mode!=null,"Ranking binding requires a recipe and explicit verification policy");
		final boolean enabled=source.containsKey("heuristic");
		for(String key:RANKING){
			require(source.containsKey(key)==enabled,"Heuristic recipe requires all ranking fields together: "+key);
		}
		final boolean gc=source.containsKey("familygc");
		for(String key:GC_RANKING){
			require(source.containsKey(key)==gc,"GC recipe requires its complete optional resource bundle: "+key);
		}
		require(!gc || enabled,"GC ranking requires the released heuristic resource bundle");
		if(gc){
			final int families=Integer.parseInt(source.get("familygcfamilies"));
			final long members=Long.parseLong(source.get("familygcmembers"));
			final double coefficient=Double.parseDouble(source.get("familygclambda"));
			require(families>0 && members>=families,"GC recipe must declare a complete nonempty family/member population");
			require(Double.isFinite(coefficient) && coefficient>=0,"GC recipe coefficient must be finite and nonnegative");
			for(String key:new String[]{"familygcsha80","familygcmapsha80"}){
				final String identity=source.get(key);
				require(identity!=null && (identity.matches("[0-9a-f]{20}") ||
					mode==FamilyCalibration.VerificationMode.STRUCTURAL && identity.equals("unmeasured")),
					"GC resource requires a valid identity under "+mode+": "+key);
			}
		}
		if(!enabled){return;}
		require("hitfourth_h020_c150".equals(source.get("heuristic")),
			"A packaged heuristic must use the released common pivot");
		final int min=Integer.parseInt(source.get("minqueryaa")),max=Integer.parseInt(source.get("maxqueryaa"));
		require(min>=1 && max>=min && max<=200,"Packaged ranking requires explicit query bounds within1..200");
		for(String key:new String[]{"lengthmedianssha80","familyhydrosha80"}){
			final String identity=source.get(key);
			require(identity!=null && (identity.matches("[0-9a-f]{20}") ||
				mode==FamilyCalibration.VerificationMode.STRUCTURAL && identity.equals("unmeasured")),
				"Packaged member metadata requires a valid declared identity under "+mode+": "+key);
		}
	}

	/** Small copied policy map supports truthful evaluator reporting of inherited package settings. */
	HashMap<String,String> rankingOptions(){
		final HashMap<String,String> result=new HashMap<String,String>();
		if(compactFamily!=null){
			for(String key:new String[]{"heuristic","minqueryaa","maxqueryaa",
					"familygclambda","gcfilter","gc_z"}){
				if(options.containsKey(key)){result.put(key,options.get(key));}
			}
			return result;
		}
		for(String key:RANKING){if(options.containsKey(key)){result.put(key,options.get(key));}}
		for(String key:GC_RANKING){if(options.containsKey(key)){result.put(key,options.get(key));}}
		return result;
	}

	/** Binds one frozen rank model to a classifier loaded from this package.
	 * This does not enable it: callers select an explicit ranking mode afterward.
	 */
	public FamilyLengthScore bindLengthRanking(FamilyTailCascade.LoadedClassifier classifier,Path model,
			String modelPin,Path medians,String medianPin) throws Exception{
		return bindLengthRanking(classifier,model,modelPin,medians,medianPin,
			FamilyCalibration.VerificationMode.SHA80,null);
	}

	/** Paired structural evaluation retains one explicit tracker through every ranking input. */
	FamilyLengthScore bindLengthRanking(FamilyTailCascade.LoadedClassifier classifier,Path model,
			String modelPin,Path medians,String medianPin,FamilyCalibration.VerificationMode mode,
			StructuralInputs lengthInputs) throws Exception{
		require(mode==verificationMode && (mode==FamilyCalibration.VerificationMode.STRUCTURAL)==(lengthInputs!=null),
			"Length-model binding mode must match its package and structural tracker");
		require(classifier!=null && options.containsKey("ddlprofiles") && options.containsKey("ddlprofilessha80"),
			"Length ranking requires this package's loaded classifier and pinned replacement profiles");
		verify();
		return new FamilyLengthScore(classifier.config(),java.nio.file.Paths.get(options.get("ddlprofiles")),
			options.get("ddlprofilessha80"),model,modelPin,medians,medianPin,mode,lengthInputs);
	}

	/** Stable resolved recipe for a reproducible run; resource controls stay outside it. */
	ByteBuilder recipe(){
		return recipe(null,null);
	}

	/** Replaces only the final gate-record component; ranking dependencies stay pinned. */
	ByteBuilder recipe(Path references,String referencePin){
		return recipe(references,referencePin,-1);
	}
	/** Records the effective shortlist without rewriting the frozen source recipe. */
	ByteBuilder recipeWithTop(int top){
		require(top>0,"Recorded shortlist cap must be positive");
		return recipe(null,null,top);
	}
	private ByteBuilder recipe(Path references,String referencePin,int top){
		return recipe(references,referencePin,top,-1);
	}
	/** Resolved recipe records the same independent fixed overrides actually supplied to the loader. */
	ByteBuilder recipeWithOverrides(int top,int ddlReject){
		require(top>=0 && ddlReject>=-1,"Recorded overrides must use the loader's cap and raw-hit floor domains");
		return recipe(null,null,top,ddlReject);
	}
	private ByteBuilder recipe(Path references,String referencePin,int top,int ddlReject){
		require(compactFamily==null,"Compact packages do not emit legacy classifier.config");
		require(references==null ? referencePin==null : Files.isRegularFile(references) && referencePin!=null &&
			(referencePin.matches("[0-9a-f]{20}") || structuralInputs!=null && referencePin.equals("unmeasured")),
			"A new reference component needs its file and sha80, or explicit unmeasured identity in structural mode");
		final ByteBuilder out=new ByteBuilder("# package=").append(path.toString()).nl();
		out.append("# package_sha80=").append(fingerprint).nl();
		if(structuralInputs!=null){out.append("verification=structural").nl();}
		for(String key:new TreeSet<String>(options.keySet())){
			if(key.equals("matrixfile") || key.equals("matrixfilesha80")){continue;}
			final String value=references!=null && key.equals("familyreferences") ?
				references.toAbsolutePath().normalize().toString() :
				references!=null && key.equals("familyreferencessha80") ? referencePin :
				top>0 && key.equals("top") ? Integer.toString(top) :
				ddlReject>=0 && key.equals("ddlreject") ? Integer.toString(ddlReject) : options.get(key);
			out.append(key).append('=').append(value).nl();
		}
		if(options.containsKey("matrix") && !options.containsKey("matrixbinding")){
			out.append("matrixbinding=").append(matrix.binding()).nl();
		}
		return out;
	}

	/** Default mode rechecks content identities. Structural mode checks observed
	 * metadata for recipe, direct components and loaded nested profiles; it cannot
	 * establish byte identity or detect changes hidden by identical metadata. */
	public void verify() throws Exception{
		if(compactFamily!=null){
			structuralInputs.verify();
			return;
		}
		if(gcInputs!=null){gcInputs.verify();}
		if(structuralInputs!=null){structuralInputs.verify();return;}
		require(fingerprint.equals(DigestSuffix.file(path.toString())),"Classifier recipe changed during evaluation");
		for(String key:PATHS){
			if(options.containsKey(key)){
				require(options.get(key+"sha80").equals(DigestSuffix.file(options.get(key))),
					"Classifier component changed during evaluation: "+key);
			}
		}
	}

	private static void require(boolean valid,String why){if(!valid){throw new IllegalArgumentException(why);}}
	final Path path;
	final String fingerprint;
	private final Path compactFamily,compactRoot,compactReference;
	private final FamilyCompactCodec.Header compactHeader;
	/** Assigned once during compact construction; never exposed without copying the array. */
	private Path[] compactInputs;
	private final ProteinSubstitutionMatrix matrix;
	private final FamilyCalibration.VerificationMode verificationMode;
	private final StructuralInputs structuralInputs;
	private final StructuralInputs gcInputs;
	/** Provenance of the fitted members; does not restrict evaluation-query lengths. */
	public final MemberLengthRange memberRange;
	private final HashMap<String,String> options=new HashMap<String,String>();
	private static final String[] PATHS={"familyreferences","ddlprofiles","scalarbounds","blacklist",
		"lengthmedians","familyhydro","matrixfile","familygc","familygcmap"};
	private static final String[] RANKING={"heuristic","lengthmedians","lengthmedianssha80","familyhydro",
		"familyhydrosha80","minqueryaa","maxqueryaa"};
	private static final String[] GC_RANKING={"familygc","familygcsha80","familygcmap","familygcmapsha80",
		"familygcfamilies","familygcmembers","familygclambda"};
	private static final String[] REQUIRED={"familyreferences","top","k","seed","width","ddlreject","ddlaccept",
		"cosine","rankingfeatures","scalargates","identityacceptfloor","blosumacceptfloor","indexearlyreject"};
	private static final String COMPACT_REFERENCE_MAGIC="#bbgap_compact_classifier_reference",
		COMPACT_REFERENCE_VERSION="1";

	private static final class CompactSource{
		CompactSource(Path family_,Path reference_){
			family=family_; reference=reference_;
		}

		final Path family,reference;
	}
}
