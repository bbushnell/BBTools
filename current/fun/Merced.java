package fun;

import java.util.HashMap;

/**
 * Demonstration class for HashMap usage and basic object-oriented programming
 * concepts. Provides educational examples of key-value data structures and
 * nested class design patterns in Java.
 * @author Brian Bushnell
 */
public class Merced {

	/**
	 * Main method demonstrating HashMap creation and type-specific mappings.
	 * Contains commented examples showing proper HashMap usage and common
	 * compilation errors.
	 * @param args Command line arguments (unused)
	 */
	public static void main(String[] args) {
		
		//Create a new Map, mapping Integers to Strings
		HashMap<String, Integer> map=new HashMap<String, Integer>();
		//Associate value 21 with key "Age"
//		map.put("Age", 21);
//		map.put("Weight", 130);
//		//This won't compile since 17.0 is a Double, not an int
//		map.put("Age", 17.0);
//		map.put("Age", 23);
//		//This won't compile since "Tall" is a String, not int
//		map.put("Height", "Tall");
//		//What will this return?
//		Integer age=map.get("Age");
//		//This will succeed, but return null ("None" in Python)
//		Integer ssn=map.get("SSN#");
//		//Won't compile because the return is an Integer
//		String weight=map.get("Weight");
	}
	
	/** Inner class representing a person with demographic attributes and
	 * computational methods for calculating derived properties. */
	class Person {
		
		/** Constructs a Person with the specified name.
		 * @param name_ The person's name */
		Person(String name_){
			name=name_;
		}
		
		/**
		 * Calculates a person's attractiveness metric based on looks, salary, and age.
		 * Uses the formula: (looks * salary) / age
		 * @return Computed attractiveness value as a float
		 */
		float attractiveness() {
			return looks*salary/age;
		}
		
		/** Person's name */
		String name;
		/** Person's age in years */
		int age;
		/** Person's salary amount */
		int salary;
		/** Person's physical appearance rating */
		float looks;
		/** Person's intelligence rating (marked for removal) */
		float brains;//TODO: Irrelevant, remove
		
	}
	
	/** HashMap mapping person names to Person objects */
	HashMap<String, Person> nameTable=new HashMap<String, Person>();
	
}
