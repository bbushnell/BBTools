package shared;

import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

/**
 * Utility to automatically copy assignable, same-named fields from a source to
 * a destination.  Designed to reduce boilerplate in constructors.  Object
 * values are copied shallowly: arrays, lists, and other references are aliased
 * rather than cloned.
 * 
 * @author Brian Bushnell
 * @contributor Sucrose
 * @date February 8, 2026
 */
public class FieldCopier{
	
	/** 
	 * Copies primitive and String fields from source to destination.  This
	 * default mode walks superclasses on both objects, skips static and final
	 * destination fields, ignores object-valued fields, and returns 0 when
	 * either input is null.
	 * @param dest Object to receive values
	 * @param source Object providing values
	 * @return Number of fields set
	 */
	public static int copy(Object dest, Object source){
		return copy(dest, source, false, false, true);
	}

	/**
	 * Copies fields from source to destination when names match and the
	 * destination type is assignable from the source type.  Reflection or access
	 * failures on one field are swallowed, allowing a partial copy.  When
	 * superclass walking is enabled, duplicate source names assert; with
	 * assertions disabled, the most-derived source declaration wins.
	 * @param dest Object to receive values
	 * @param source Object providing values
	 * @param includeObjects If true, copies complex Objects (Arrays, Lists, etc). 
	 * If false, only copies Primitives and Strings.
	 * @param includeStatics If true, also attempts to copy static fields.
	 * @param includeSuperclass If true, walks parent classes on both objects.
	 * @return Number of fields set
	 */
	public static int copy(Object dest, Object source, 
			boolean includeObjects, boolean includeStatics, boolean includeSuperclass){
		if(dest==null || source==null){return 0;}
		
		// Map destination fields for O(1) lookup.
		Map<String, Field> destFields=getAllFields(dest.getClass(), includeSuperclass);
		Class<?> sourceClass=source.getClass();
		int matches=0;
		final HashSet<String> seenSourceNames=new HashSet<String>();//[FieldCopier#001] names already claimed by a more-derived declaration

		while(sourceClass!=null && sourceClass!=Object.class){
			final Field[] fields=sourceClass.getDeclaredFields();
			
			for(Field srcField : fields){
				try{
					final String name=srcField.getName();
					final boolean fresh=seenSourceNames.add(name);
					assert(fresh) : "FieldCopier: shadowed field '"+name+"' in "+source.getClass()+
						" (same name declared at multiple hierarchy levels).";
					if(!fresh){continue;}
					final int srcMods=srcField.getModifiers();
					
					if(Modifier.isStatic(srcMods) && !includeStatics){continue;}
					
					final Class<?> type=srcField.getType();
					final boolean isPrimOrString=(type.isPrimitive() || type==String.class);
					if(!isPrimOrString && !includeObjects){continue;}
					
					final Field destField=destFields.get(name);
					if(destField==null){continue;}
					
					final int destMods=destField.getModifiers();
					if(Modifier.isFinal(destMods)){continue;}
					if(Modifier.isStatic(destMods) && !includeStatics){continue;}

					if(destField.getType().isAssignableFrom(type)){
						if(!srcField.isAccessible()){srcField.setAccessible(true);}
						if(!destField.isAccessible()){destField.setAccessible(true);}
						
						final Object val=srcField.get(source);
						destField.set(dest, val);
						matches++;
					}
				}catch(Exception e){
					// Security exceptions or access violations are ignored.
				}
			}
			sourceClass=(includeSuperclass ? sourceClass.getSuperclass() : null);
		}
		return matches;
	}
	
	/**
	 * Returns declared fields in a class, optionally walking up the hierarchy.
	 * Handles shadowing by prioritizing the lowest class in the hierarchy.
	 * @param clazz Initial class to inspect
	 * @param recursive True to include ancestor classes
	 * @return Field map by name
	 */
	private static Map<String, Field> getAllFields(Class<?> clazz, boolean recursive){
		Map<String, Field> map=new HashMap<>();
		while(clazz!=null && clazz!=Object.class){
			for(Field f : clazz.getDeclaredFields()){
				map.putIfAbsent(f.getName(), f);
			}
			clazz=(recursive ? clazz.getSuperclass() : null);
		}
		return map;
	}
}
