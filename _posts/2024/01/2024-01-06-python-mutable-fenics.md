---
layout: post
title: "Understanding Object References in Python: Immutable vs Mutable Behaviors"
description: Exploring Mutable and Immutable Object References in Python and Their Implications in FEniCS Simulations.
categories:
  - coding
tags:
  - fenics
  - python
---
In Python, the way variables are handled is based on what they are: objects. Everything in Python is an object, including integers, floats, functions, and even types and classes. However, the behavior of these objects can vary significantly depending on whether they are mutable or immutable.

## Immutable Objects:
Immutable objects in Python include `integers`, `floats`, `strings`, and `tuples`. 
- When you create an immutable object, like setting `u = 10`, you create an object `10`, and `u` is a name that refers to it. If you then do `v = u`, you are creating a new reference `v` that also points to object `10`. 
- However, because integers are immutable when you do something like `u += 1`, you are not modifying the object `10` — you are creating a new object `11` and changing `u` to refer to it instead. This is a critical distinction: `u` and `v` are separate references, and changing one does not affect the other.

```python
u = 10
v = u  # v is now a reference to the same integer object as u

print("Original u:", u)  # Output: 10
print("Original v:", v)  # Output: 10

u += 1  # u is now a reference to a new integer object, 11
print("Modified u:", u)  # Output: 11
print("Unchanged v:", v)  # Output: 10

```

When you modify `u` by doing `u += 1`, you're not changing the integer object `10`. Instead, you're creating a new integer object `11` and updating `u` to reference this new object. Since integers are immutable, this is the only way to "change" `u`. The variable `v` still references the original integer object `10`.
## Mutable Objects:

Mutable objects include `lists`, `dictionaries`, and most `class instances`. The word "most" acknowledges the existence of deliberately designed immutable classes in Python, even though they are less common. Most classes you'll encounter or write are likely mutable because that's the default behavior and often what's desired. However, when designing or using a class, it's important to understand its mutability, as it will affect how you interact with its instances and what implications it has for your program's behavior.

With mutable objects such as instances of custom classes, assignments don't create new objects; they only create new references to the same object.

```python
class CustomObject:
    def __init__(self, value):
        self.value = value

u = CustomObject(10)
v = u  # v is now another reference to the same CustomObject instance as u

print("Original u.value:", u.value)  # Output: 10
print("Original v.value:", v.value)  # Output: 10

u.value += 1  # Modifies the value attribute of the object referred to by u
print("Modified u.value:", u.value)  # Output: 11
print("v.value also changes:", v.value)  # Output: 11

```

In this mutable scenario, when you modify `u.value` by doing `u.value += 1`, you modify the object that both `u` and `v` reference. As a result, the change is reflected in both `u.value` and `v.value`.

## A note on FEniCS
Suppose you have a FEniCS `Function` object `u` and assign it to an attribute of a class instance like `self.u`. In this case, both `u` and `self.u` refer to the same actual object in memory. If this object is mutable, any changes made to it through any of its references will be reflected across all references. For instance, if the function object `u` is modified in place in a method of your class, this change will also be visible outside of the class because `self.u` and `u` are just two names for the same object.

This also holds true for passing a FEniCS `Function` object to a python function.
```python
from fenics import *

# Define the problem domain, mesh, function space
mesh = UnitSquareMesh(8, 8)
V = FunctionSpace(mesh, 'P', 1)

# Create a FEniCS Function object (mutable)
u = Function(V)

def modify_function(fn):
    # Modify the function object
    fn.vector()[:] += 1  # Increment all values by 1
    print("Values of function inside modify_function:", fn.vector()[:])

print("Values of u before function call:", u.vector()[:])
modify_function(u)  # Pass the mutable object to the function
print("Values of u after function call:", u.vector()[:])  # u is modified outside function

```

In this FEniCS example, `u` is a `Function` object that is mutable. When passed to the `modify_function`, the function modifies the `Function` object `fn` (the same as `u`). The modification is an in-place change, incrementing all function values by 1. After the function call, you can observe that `u` has been modified outside the function because `fn` and `u` refer to the same `Function` object.

### Further discussion on the assignment in FEniCS

1. `u.assign(v)`: This method is specifically used in FEniCS. When you have two function objects, say `u` and `v`, `u.assign(v)` assigns the values of `v` to `u` without changing the underlying function space of `u`. This means that after the assignment, `u` and `v` will have the same values in their respective function spaces, but their function spaces or any other associated data remain unchanged. This method is particularly used when you want to update a solution or a variable iteratively without altering its properties.

2. `u = v`: This is a standard assignment operation in Python. In the context of FEniCS, when you say `u = v`, you are not just copying the values from `v` to `u`. Still, you are also making `u` reference the same object as `v`. After this operation, any change to `v` will reflect in `u` as well because they are essentially the same object. This doesn't just assign the current values but makes both variables refer to the same underlying data and function space.

In summary, `u.assign(v)` is the way to copy values from one function to another in FEniCS, preserving the original structure and properties of the function, while `u = v` is a way to make two variables reference the same object, leading to a complete overlap of their identities. 

In practice, `u.assign(v)` is used when you want to update or iterate solutions in FEniCS, keeping the function spaces distinct, whereas `u = v` is a broader Python assignment operation, which in FEniCS context would typically be used for setting up initial conditions or simplifying references.

## Assignment Summary:
- **Immutable objects**: Assignment creates a new reference to an object. If you modify the original reference (like changing `u` to point to a different integer), the new reference (`v`) is unaffected because it still points to the original object.
- **Mutable objects**: Assignment creates a new reference to the same object. If you modify the object through one reference, all references see the change because they all point to the same object.

Understanding these behaviors is crucial when working with different types of objects in Python, especially when the effects of functions and methods on objects can significantly impact the program's behavior. It's also why careful management of object references is important in Python to avoid unintended side effects, particularly with mutable objects.
