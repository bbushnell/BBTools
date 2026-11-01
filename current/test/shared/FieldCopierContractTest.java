package test.shared;

import shared.FieldCopier;

/** Finite FieldCopier contract checks for static/final filters and hierarchy handling.
 * Cases use literal independent expectations and no production callers.
 * @author Jean
 */
public final class FieldCopierContractTest{

	/** Runs every case and fails after reporting all case results. */
	public static void main(String[] args){
		if(args.length!=0){throw new IllegalArgumentException("Expected no arguments");}
		final FieldCopierContractTest test=new FieldCopierContractTest();
		test.check("default-static-destination-skip", FieldCopierContractTest::defaultSkipsStaticDestinations);
		test.check("include-object-control", FieldCopierContractTest::includeObjectsCopiesAssignableReference);
		test.check("default-superclass-control", FieldCopierContractTest::defaultIncludesSuperclasses);
		test.check("exclude-superclass-control", FieldCopierContractTest::excludeSuperclassesCopiesOnlyLeaf);
		test.check("exclude-source-superclass", FieldCopierContractTest::excludeSourceSuperclassFields);
		test.check("include-static-control", FieldCopierContractTest::includeStaticsCopiesStaticFields);
		System.out.println("FieldCopierContractTest cases="+test.cases+" failures="+test.failures);
		if(test.failures!=0){throw new AssertionError("FieldCopier contract failures="+test.failures);}
	}

	/** Default copy must skip static and final destinations, while copying primitives and Strings. */
	private static void defaultSkipsStaticDestinations(){
		StaticDest.staticNumber=3;
		final StaticDest dest=new StaticDest();
		final StaticSource source=new StaticSource();
		final int copied=FieldCopier.copy(dest, source);
		require(copied==2, "default copied "+copied+" fields instead of instance+text only");
		require(StaticDest.staticNumber==3, "default copy changed a static destination");
		require(dest.instanceNumber==11, "instance primitive was not copied");
		require("copied".equals(dest.text), "String was not copied");
		require(dest.payload==null, "Object payload was copied by default");
		require(dest.finalNumber==5, "final destination changed");
	}

	/** Object fields are copied only when explicitly requested and assignable. */
	private static void includeObjectsCopiesAssignableReference(){
		final ObjectDest dest=new ObjectDest();
		final ObjectSource source=new ObjectSource();
		final int copied=FieldCopier.copy(dest, source, true, false, true);
		require(copied==1, "includeObjects copied "+copied+" fields instead of one");
		require(dest.payload==source.payload, "assignable object reference was not copied");
	}

	/** The default overload walks parent fields on both source and destination. */
	private static void defaultIncludesSuperclasses(){
		final ChildDest dest=new ChildDest();
		final ChildSource source=new ChildSource();
		final int copied=FieldCopier.copy(dest, source);
		require(copied==2, "default superclass copy count="+copied);
		require(dest.parent==21 && dest.child==34, "superclass or child primitive was not copied");
	}

	/** Disabling superclass traversal restricts copying to the leaf class. */
	private static void excludeSuperclassesCopiesOnlyLeaf(){
		final ChildDest dest=new ChildDest();
		final ChildSource source=new ChildSource();
		final int copied=FieldCopier.copy(dest, source, false, false, false);
		require(copied==1, "leaf-only copy count="+copied);
		require(dest.parent==0 && dest.child==34, "leaf-only copy touched a parent or missed the child");
	}

	/** Disabling superclass traversal must also skip fields inherited only by the source. */
	private static void excludeSourceSuperclassFields(){
		final LeafOnlyDest dest=new LeafOnlyDest();
		final SourceChildWithOnlyParent source=new SourceChildWithOnlyParent();
		final int copied=FieldCopier.copy(dest, source, false, false, false);
		require(copied==0, "leaf-only copy used a source parent field");
		require(dest.parentOnly==0, "leaf-only copy changed destination from a source parent");
	}

	/** Explicit static inclusion copies matching static source fields to static destinations. */
	private static void includeStaticsCopiesStaticFields(){
		StaticOnlyDest.value=7;
		StaticOnlySource.value=99;
		final int copied=FieldCopier.copy(new StaticOnlyDest(), new StaticOnlySource(), false, true, true);
		require(copied==1, "includeStatics copy count="+copied);
		require(StaticOnlyDest.value==99, "explicit static field was not copied");
	}

	/** Runs one case and records rather than short-circuiting failures. */
	private void check(String name, Action action){
		cases++;
		try{
			action.run();
			System.out.println(name+" PASS");
		}catch(AssertionError e){
			failures++;
			System.out.println(name+" FAIL: "+e.getMessage());
		}
	}

	/** Checks expectations in both assertion modes. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** One finite FieldCopier contract case. */
	private interface Action{void run();}

	/** Destination with mutable static, mutable instance, final, String, and Object fields. */
	private static class StaticDest{
		static int staticNumber=3;
		int instanceNumber=0;
		final int finalNumber=5;
		String text="original";
		Object payload=null;
	}

	/** Source with matching instance fields; staticNumber must not hit StaticDest.staticNumber by default. */
	private static class StaticSource{
		int staticNumber=17;
		int instanceNumber=11;
		int finalNumber=19;
		String text="copied";
		StringBuilder payload=new StringBuilder("object");
	}

	/** Object reference destination. */
	private static class ObjectDest{Object payload;}
	/** Assignable object reference source. */
	private static class ObjectSource{StringBuilder payload=new StringBuilder("payload");}

	/** Destination parent field. */
	private static class ParentDest{int parent;}
	/** Destination leaf field. */
	private static class ChildDest extends ParentDest{int child;}
	/** Source parent field. */
	private static class ParentSource{int parent=21;}
	/** Source leaf field. */
	private static class ChildSource extends ParentSource{int child=34;}

	/** Leaf destination with a field that exists only in the source parent. */
	private static class LeafOnlyDest{int parentOnly;}
	/** Source parent with the only matching field. */
	private static class SourceParentOnly{int parentOnly=55;}
	/** Source leaf deliberately has no matching fields. */
	private static class SourceChildWithOnlyParent extends SourceParentOnly{}

	/** Static-only destination. */
	private static class StaticOnlyDest{static int value=7;}
	/** Static-only source. */
	private static class StaticOnlySource{static int value=99;}

	/** Attempted cases. */
	private int cases;
	/** Failed cases. */
	private int failures;
}
