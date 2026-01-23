# Contributing to DREAM

Thank you for your interest in contributing to DREAM!

DREAM is a long-lived scientific software project developed and maintained by a distributed community of researchers, students, and developers. Over the years, it has grown to support a wide range of physics models and use cases, and contributions - large and small - have played an important role in that evolution.

Contributions to DREAM can take many forms, including:
- Bug fixes and robustness improvements
- New physics models or numerical methods
- Performance optimizations
- Improvements to documentation, examples, or tests
- Refactoring or cleanup that improves clarity or maintainability

Both new contributors and returning contributors are very welcome. If you are unsure where to start, feel free to open an issue to ask questions or discuss ideas before writing code.

---

## General contribution philosophy

DREAM is primarily a **scientific computing codebase**, and contributions are evaluated with an emphasis on:
- Physical correctness and numerical robustness
- Clarity and maintainability
- Performance awareness in computationally intensive parts of the code

The codebase has evolved over many years, and not all parts follow a single uniform style. Rather than enforcing strict rules, the goal of this guide is to document conventions that have worked well and to help contributors write code that fits naturally into the existing project.

In particular:
- Small, focused pull requests are preferred over large, sweeping changes.
- Following the local style of the code you are modifying is usually the best choice.
- If a contribution involves significant design decisions or API changes, it is encouraged to discuss them in an issue first.

The sections below provide practical guidance intended to make contributing smoother and reviews more predictable.

---

## Coding style and guidelines


This project has grown over the years and contains a mix of styles. The goal of these guidelines is **not** to enforce uniformity everywhere, but to document the conventions that have generally worked well and to help contributors write code that fits naturally into the existing codebase.

When contributing new code or modifying existing code, the most important rule is:

> **Follow the style of the surrounding code.**

These guidelines are meant to provide orientation and shared expectations, not to impose rigid rules.

---

### What this guide is not

To set expectations clearly, it may be helpful to state what these guidelines are **not** intended to do:

- **This is not a strict rulebook.**  
  The guidelines are descriptive and permissive. They document conventions that have generally worked well in the project, rather than enforcing a single “correct” way to write code.

- **This is not an attempt to retroactively standardize the entire codebase.**  
  Existing code is not expected to be rewritten solely to conform to these guidelines. Incremental improvements are welcome, but large mechanical refactors should be motivated by clear benefits beyond stylistic consistency.

- **This is not a replacement for good judgment.**  
  Scientific and performance-sensitive code often involves trade-offs. In some cases, deviating from these guidelines may be the right choice; such deviations are acceptable when they improve clarity, correctness, or performance.

- **This is not meant to introduce unnecessary barriers to contribution.**  
  The primary goal is to make contributions easier, reviews smoother, and expectations clearer - not to add process or friction.

When in doubt, contributors are encouraged to prioritize clarity, follow the local style of the code they are working in, and ask for feedback during review.

---

### General principles

- Prefer **clarity over cleverness**, especially in numerically sensitive code.
- Make **data flow, memory layout, and ownership explicit**.
- Favor **simple, explicit C++** over advanced language features unless there is a clear benefit.
- Be mindful of performance and allocation behavior in hot paths.
- When in doubt, consistency with existing code is more important than strict adherence to any single rule.

---

### C++ language usage

- The codebase intentionally uses a **C-leaning subset of C++**.
- Avoid heavy template metaprogramming, complex type traits, or implicit abstractions unless there is a strong justification.
- Use STL containers (`std::vector`, `std::map`, etc.) where appropriate, but be mindful of their cost in performance-critical sections.
- Prefer explicit loops and indexing over algorithm abstractions in hot paths.
- Use modern C++ features sparingly and only when they improve clarity or safety without obscuring intent.

---

### Memory management

- Manual memory management (`new[]`, `delete[]`) is used in several parts of the codebase and is acceptable when:
  - Ownership is clear
  - Allocation and deallocation are symmetric
  - Object lifetime is well-defined
- Functions or classes that take ownership of memory should **clearly document this behavior**.
- When memory is shared or externally owned, this should be made explicit (e.g. through naming, flags, or documentation).
- Avoid hidden or implicit allocations in performance-critical code.

---

### Naming conventions

Naming conventions are intentionally flexible and follow existing patterns:

- Classes and public APIs generally use `CamelCase`.
- Member functions typically use `CamelCase`.
- Local variables use descriptive `lowerCamelCase` names, or short names (`i`, `j`, `nr`) where the domain meaning is clear.
- Constants and enum values use `ALL_CAPS`.
- Type aliases follow established conventions (e.g. `real_t`, `len_t`, `int_t`).

Consistency with surrounding code is more important than strict adherence to any single naming rule.

---

### File and code organization

Source files typically follow this structure:

1. File-level comment
2. Includes
3. Namespace usage
4. Constants and static definitions
5. Constructors and destructors
6. Public API methods
7. Internal helpers

- Long functions are acceptable when they represent a single, coherent algorithm.
- Logical sections within functions should be separated by whitespace and/or comments.

---

### Comments and documentation

- Use comments to explain **assumptions, limitations, or non-obvious behavior**.
- Comments should focus on *why* something is done, not restating the code.
- Functions and classes should be documented using block comments (`/** ... */`) where appropriate.
- Inline comments are encouraged when code relies on domain knowledge or makes non-obvious decisions.
- Explicit markers such as `TODO` or `XXX` are acceptable to flag known limitations or future work.

---

### Error handling

- Prefer explicit runtime checks and informative error messages.
- Fail fast when encountering invalid states or inconsistent input.
- Error messages should provide enough context to diagnose scientific or numerical issues.

---

## Tests and documentation

Contributions are encouraged to include appropriate tests and documentation updates where relevant. The level and type of testing will naturally depend on the nature of the change.

In particular:
- **Bug fixes** should ideally include a regression test that would have failed before the fix.
- **New functionality or models** are encouraged to come with tests that validate expected behavior.
- **Numerical or physics-related changes** may be best validated through integration or physics-level tests rather than isolated unit tests.

DREAM contains a mix of:
- Unit tests for individual components
- Integration or “physics” tests that validate coupled behavior or physical assumptions

Contributors are encouraged to use the testing approach that best matches the change being made.

For user-facing changes, additions to the **Sphinx documentation**, docstrings, or examples are encouraged so that new functionality is discoverable and usable by others.

If you are unsure what level of testing or documentation is appropriate, feel free to ask in an issue or pull request - reviewers are happy to help clarify expectations.
