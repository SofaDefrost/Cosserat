1. Prepare a clean working state:

   - Save or stash any uncommitted changes in the current working directory.
   - Pull the latest changes from remote_master to ensure you have the most up-to-date code base.

2. Create and switch to the feature branch:

   - Create a new feature branch named "feature/lie-groups" from remote_master.
   - Confirm that you are on the new branch before proceeding with any code changes.

3. Set up the directory structure:

   - Inside the Cosserat directory, create a subdirectory named "liegroups/".
   - Add the following initial files:
     • liegroups/LieGroupBase.h (templated base class declaration)
     • liegroups/LieGroupBase.inl (template implementations)
     • liegroups/Types.h (type definitions and forward declarations)

4. Implement the LieGroupBase class:

   - Declare a templated class LieGroupBase to define common methods for Lie groups.
   - Provide pure virtual (or abstract) functions for composition, inverse, exponential, logarithm, Adjoint, and group action.
   - Use Eigen for matrix/vector operations.

5. Implement specific derived classes:

   - Create derived classes for ℝ(n), SO(2), SE(2), SO(3), SE(3), SE_2(3), and SGal(3) by extending LieGroupBase.
   - Implement each class’s specific group properties and overrides of the base class functions.
   - Leverage Eigen’s transforms (e.g., Rotation2D, Quaternion, Isometry3d) where appropriate.

6. Implement the Bundle template class:

   - Add a Bundle class to liegroups/ for handling manifold products or bundles.
   - Provide functionalities to combine multiple Lie groups into a unified manifold representation.

7. Integrate into the build and project structure:

   - Update CMakeLists.txt to include the new liegroups/ directory and files in the build process.
   - Ensure the created headers are referenced correctly in Cosserat/fwd.h and any other relevant headers.
   - Follow the existing code style and documentation guidelines throughout.

8. Test and verify:

   - Compile and run tests (existing or newly added) to validate the new Lie group functionalities.
   - Review code style compliance and documentation for clarity and consistency.

9. Finalize and push:
   - Commit the changes with an informative commit message.
   - Push the new feature branch to the remote repository for review or pull request.
